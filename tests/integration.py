#!/usr/bin/env python3
"""End-to-end importer checks in disposable databases on an explicit test server.

python3 tests/integration.py --binary parser/build/cif_parser --dsn 'host=... dbname=postgres'
Never uses or alters the application's unncdb database.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import uuid

ROOT = Path(__file__).resolve().parents[1]


def execute(dsn, sql):
    result = subprocess.run(['psql', '-X', '-qAt', '-v', 'ON_ERROR_STOP=1', dsn],
                            input=sql, text=True, capture_output=True, check=True)
    return result.stdout.strip()


def schema(dsn, path):
    subprocess.run(['psql', '-X', '-q', '-v', 'ON_ERROR_STOP=1', dsn, '-f', str(path)],
                   check=True, capture_output=True, text=True)


def run_import(binary, path, dsn, *flags, expected=0):
    r = subprocess.run([str(binary), str(path), dsn, *flags], capture_output=True, text=True)
    if r.returncode != expected:
        raise AssertionError(f'import exit {r.returncode}, wanted {expected}: {r.stderr}')
    return json.loads(r.stdout) if r.returncode == 0 else r


def block(name, atoms='C1 C 0.1 0.2 0.3', extra='', occupancy=False):
    return f'''data_{name}
_CELL_LENGTH_A 10
_cell_length_b 11
_cell_length_c 12
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_chemical_formula_moiety 'C1'
_symmetry_Int_Tables_number 2
_space_group_name_Hall '-P 1'
loop_
_symmetry_equiv_pos_as_xyz
'x,y,z'
'-x,-y,-z'
{extra}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{'_atom_site_occupancy' if occupancy else ''}
{atoms}
'''.encode()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', required=True, type=Path)
    ap.add_argument('--dsn', required=True, help='Maintenance connection to a TEST PostgreSQL server')
    args = ap.parse_args()
    binary = args.binary.resolve()
    names = ['unncdb_test_' + uuid.uuid4().hex]
    created = []
    try:
        for name in names:
            execute(args.dsn, f'CREATE DATABASE {name}')
            created.append(name)
        fresh = args.dsn + ' dbname=' + names[0]
        schema(fresh, ROOT/'unncdb_schema.sql')
        with tempfile.TemporaryDirectory(prefix='unncdb-integration-') as td:
            td = Path(td)
            path = td/'mixed.cif'
            source = b'# exact preamble\r\n' + block('good', extra="_refine_special_details\n;\ndata_fake_inside_text\n;\n")
            source += block('problem', 'C1 C 1.23(4)e-2 0 0 1\nC1 N ? 0 0 ?\nCa1? Ca 0 0 0 .',
                            extra="loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\nC1 0.02\n", occupancy=True)
            source += b'data_bad\n_tag_without_value\n'
            source += b'data_binary\n_note "bad\x00\xff"\n'
            source = source.replace(b'\n', b'\r\n')  # Preserve even unusual CR/CRLF combinations exactly.
            path.write_bytes(source)
            first = run_import(binary,path,fresh,'--batch-size','2')
            assert first['blocks_seen'] == first['blocks_inserted'] == 4, first
            assert first['atom_rows_extracted'] == 4, first
            assert execute(fresh, 'SELECT count(*) FROM atoms') == '4'
            assert execute(fresh, 'SELECT count(*) FROM atom_adp WHERE atom_id IS NULL') == '1'
            assert execute(fresh, "SELECT count(*) FROM atoms WHERE label='C1'") == '3'
            assert execute(fresh, "SELECT occupancy_state FROM atoms WHERE source_row=1 ORDER BY atom_id LIMIT 1") == 'absent'
            assert execute(fresh, "SELECT occupancy FROM atoms WHERE occupancy_state='absent'") == '1'
            assert execute(fresh, "SELECT count(*) FROM atoms WHERE occupancy IS NULL") == '2'
            assert execute(fresh, "SELECT count(*) FROM structures WHERE is_centrosymmetric") == '2'
            assert execute(fresh, "SELECT count(*) FROM cif_blocks WHERE parse_status='parse_error'") == '2'
            assert execute(fresh, "SELECT count(*) FROM import_issues WHERE code='duplicate_atom_label'") == '1'
            assert execute(fresh, "SELECT count(*) FROM atoms WHERE x IS NULL") == '1'
            preamble = bytes.fromhex(execute(fresh,"SELECT encode(preamble,'hex') FROM cif_sources"))
            rows = execute(fresh,"SELECT byte_offset||'|'||byte_length||'|'||encode(raw_cif,'hex')||'|'||sha256 FROM cif_blocks ORDER BY ordinal")
            restored = preamble
            for row in rows.splitlines():
                offset,length,raw,digest = row.split('|')
                raw = bytes.fromhex(raw)
                assert source[int(offset):int(offset)+int(length)] == raw
                assert hashlib.sha256(raw).hexdigest() == digest
                restored += raw
            assert restored == source
            assert execute(fresh,'SELECT sha256 FROM cif_sources') == hashlib.sha256(source).hexdigest()
            renamed = td/'same-content.cif';shutil.copyfile(path,renamed)
            again = run_import(binary,renamed,fresh)
            assert again['blocks_inserted'] == 0 and again['blocks_reused'] == 4, again
            assert execute(fresh,'SELECT count(*) FROM atoms') == '4'
            # A database error in batch two must roll it back, retaining batch one.
            rollback = td/'rollback.cif'
            rollback.write_bytes(b''.join(block(f'rollback{i}') for i in range(6)))
            execute(fresh,"""CREATE FUNCTION reject_test_block() RETURNS trigger LANGUAGE plpgsql AS $$
            BEGIN IF NEW.block_name='rollback3' THEN RAISE EXCEPTION 'injected test failure'; END IF; RETURN NEW; END $$;
            CREATE TRIGGER reject_test_block BEFORE INSERT ON cif_blocks FOR EACH ROW EXECUTE FUNCTION reject_test_block();""")
            run_import(binary,rollback,fresh,'--batch-size','2',expected=1)
            assert execute(fresh,"SELECT count(*) FROM cif_blocks WHERE block_name LIKE 'rollback%'") == '2'
            assert execute(fresh,"SELECT count(*) FROM structures WHERE refcode LIKE 'rollback%'") == '2'
            execute(fresh,'DROP TRIGGER reject_test_block ON cif_blocks; DROP FUNCTION reject_test_block()')
            resumed = run_import(binary,rollback,fresh,'--batch-size','2')
            assert resumed['blocks_reused'] == 2 and resumed['blocks_inserted'] == 4, resumed
            assert execute(fresh,"SELECT count(*) FROM structures WHERE refcode LIKE 'rollback%'") == '6'
            # Hard interruption: next run acquires the source lock and marks stale run failed.
            crash = td/'interrupted.cif'
            crash.write_bytes(b''.join(block(f'crash{i}') for i in range(1200)))
            process = subprocess.Popen([str(binary),str(crash),fresh,'--batch-size','1'],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
            assert 'Blocks:' in process.stderr.readline()
            process.kill();process.wait();process.stderr.close()
            partial = int(execute(fresh,"SELECT count(*) FROM cif_blocks WHERE block_name LIKE 'crash%'") or 0)
            assert 0 < partial < 1200, partial
            resumed = run_import(binary,crash,fresh,'--batch-size','100')
            assert resumed['blocks_reused'] == partial and resumed['blocks_seen'] == 1200, resumed
            assert execute(fresh,"SELECT count(*) FROM structures WHERE refcode LIKE 'crash%'") == '1200'
            assert execute(fresh,"SELECT count(*) FROM cif_imports WHERE status='running'") == '0'
        # Schema initialization is for empty databases: a second application
        # must fail atomically without altering already imported records.
        before = execute(fresh, 'SELECT count(*) FROM structures')
        try:
            schema(fresh, ROOT/'unncdb_schema.sql')
        except subprocess.CalledProcessError:
            pass
        else:
            raise AssertionError('Applying the schema to a populated database must fail')
        assert execute(fresh, 'SELECT count(*) FROM structures') == before
        assert execute(fresh, 'SELECT version FROM schema_migrations') == '2'
        assert execute(fresh, "SELECT pg_get_constraintdef(oid) FROM pg_constraint WHERE conrelid='atoms'::regclass AND contype='p'") == 'PRIMARY KEY (atom_id)'
        # NULL labels still represent stored atom rows and must be counted.
        execute(fresh, "INSERT INTO atoms(structure_id,source_row) SELECT id,999 FROM structures WHERE refcode='good'")
        assert execute(fresh, "SELECT n_atoms FROM structures_summary WHERE refcode='good'") == '2'
        # Both linked and ambiguous ADP rows retain proper cascading references.
        execute(fresh, """INSERT INTO atom_adp(structure_id,atom_id,source_row,label,u11)
        SELECT structure_id,atom_id,999,label,0.02 FROM atoms
        WHERE label='C1' AND structure_id=(SELECT id FROM structures WHERE refcode='good');
        DELETE FROM atoms WHERE label='C1' AND structure_id=(SELECT id FROM structures WHERE refcode='good');""")
        assert execute(fresh, "SELECT count(*) FROM atom_adp WHERE source_row=999") == '0'
        assert execute(fresh, "SELECT n_atoms FROM structures_summary WHERE refcode='good'") == '1'
        execute(fresh, "DELETE FROM structures WHERE refcode='problem'")
        assert execute(fresh, 'SELECT count(*) FROM atom_adp') == '0'
        print('PASS: standalone schema, exact bytes, invalid blocks, all atom rows, occupancy states, symmetry, retries, rollback, interrupted resume, atom identity, summary counts, ADP cascades, atomic initialization')
    finally:
        for name in reversed(created):
            execute(args.dsn,f'DROP DATABASE {name} WITH (FORCE)')


if __name__ == '__main__':
    main()
