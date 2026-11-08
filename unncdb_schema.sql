-- UNNCDB: complete schema for an empty PostgreSQL database.
-- Run: psql -X -v ON_ERROR_STOP=1 -d unncdb -f unncdb_schema.sql
-- All objects are created in one transaction. Existing tables are not replaced.
BEGIN;

-- Compatibility marker read by the C++ importer and Python module.
-- Kept at version 2 so existing installations and new databases use the same API.
CREATE TABLE schema_migrations (
    version     INTEGER PRIMARY KEY,
    applied_at  TIMESTAMPTZ NOT NULL DEFAULT now()
);

-- Original files and import runs.
CREATE TABLE cif_sources (
    id                BIGSERIAL PRIMARY KEY,
    sha256            TEXT NOT NULL UNIQUE CHECK (sha256 ~ '^[0-9a-f]{64}$'),
    byte_size         BIGINT NOT NULL CHECK (byte_size >= 0),
    original_path     TEXT NOT NULL,
    integrity_status  TEXT NOT NULL DEFAULT 'pending'
        CHECK (integrity_status IN ('pending', 'verified', 'mismatch')),
    preamble          BYTEA NOT NULL DEFAULT ''::bytea,
    created_at        TIMESTAMPTZ NOT NULL DEFAULT now()
);

CREATE TABLE cif_imports (
    id                BIGSERIAL PRIMARY KEY,
    source_id         BIGINT NOT NULL REFERENCES cif_sources(id),
    input_path        TEXT NOT NULL,
    importer_version  TEXT NOT NULL,
    parser_version    TEXT NOT NULL,
    status            TEXT NOT NULL
        CHECK (status IN ('running', 'complete', 'complete_with_issues', 'failed')),
    started_at        TIMESTAMPTZ NOT NULL DEFAULT now(),
    finished_at       TIMESTAMPTZ,
    blocks_seen       BIGINT NOT NULL DEFAULT 0,
    blocks_inserted   BIGINT NOT NULL DEFAULT 0,
    blocks_reused     BIGINT NOT NULL DEFAULT 0,
    error_message     TEXT
);

-- Preserve every block, including blocks that could not be parsed.
CREATE TABLE cif_blocks (
    id            BIGSERIAL PRIMARY KEY,
    source_id     BIGINT NOT NULL REFERENCES cif_sources(id),
    import_id     BIGINT NOT NULL REFERENCES cif_imports(id),
    ordinal       BIGINT NOT NULL CHECK (ordinal > 0),
    block_name    TEXT NOT NULL,
    byte_offset   BIGINT NOT NULL CHECK (byte_offset >= 0),
    byte_length   BIGINT NOT NULL CHECK (byte_length > 0),
    sha256        TEXT NOT NULL CHECK (sha256 ~ '^[0-9a-f]{64}$'),
    raw_cif       BYTEA NOT NULL CHECK (octet_length(raw_cif) = byte_length),
    parse_status  TEXT NOT NULL
        CHECK (parse_status IN ('parsed', 'parse_error', 'extraction_error')),
    UNIQUE (source_id, ordinal)
);

CREATE TABLE structures (
    id                        BIGSERIAL PRIMARY KEY,
    source_block_id           BIGINT REFERENCES cif_blocks(id),
    refcode                   TEXT,
    formula                   TEXT,
    formula_moiety            TEXT,
    description               TEXT,
    refinement_details        TEXT,

    -- Source symmetry descriptors and the resolved group/label.
    space_group               TEXT,
    hall_symbol               TEXT,
    space_group_num_raw       TEXT,
    declared_space_group_num  INTEGER,
    space_group_num           INTEGER,
    point_group               TEXT,
    crystal_system            TEXT,
    is_centrosymmetric        BOOLEAN,
    symmetry_status           TEXT,
    symmetry_method           TEXT,
    symmetry_version          TEXT,

    -- Unit cell: lengths in angstroms, angles in degrees.
    a                         DOUBLE PRECISION,
    b                         DOUBLE PRECISION,
    c                         DOUBLE PRECISION,
    alpha                     DOUBLE PRECISION,
    beta                      DOUBLE PRECISION,
    gamma                     DOUBLE PRECISION,
    a_su                      DOUBLE PRECISION,
    b_su                      DOUBLE PRECISION,
    c_su                      DOUBLE PRECISION,
    alpha_su                  DOUBLE PRECISION,
    beta_su                   DOUBLE PRECISION,
    gamma_su                  DOUBLE PRECISION,
    volume                    DOUBLE PRECISION,
    reduced_volume            DOUBLE PRECISION,
    z                         INTEGER,

    measurement_temperature   DOUBLE PRECISION,
    radiation_probe           TEXT,
    radiation_wavelength      DOUBLE PRECISION,
    r_factor_gt               DOUBLE PRECISION,
    r_factor_all              DOUBLE PRECISION,
    wr_factor_gt              DOUBLE PRECISION,
    wr_factor_ref             DOUBLE PRECISION,
    disorder_status           TEXT NOT NULL DEFAULT 'unknown',
    extraction_status         TEXT NOT NULL DEFAULT 'legacy',
    atom_rows_source          INTEGER,
    atom_rows_stored          INTEGER,
    source_file               TEXT,
    created_at                TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE symmetry_operations (
    structure_id   BIGINT NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    source_row     INTEGER NOT NULL,
    operation_xyz  TEXT NOT NULL,
    PRIMARY KEY (structure_id, source_row)
);

-- Atom identity is independent of its label. Invalid/missing values stay NULL;
-- the source row and import diagnostics retain the information needed to inspect them.
CREATE TABLE atoms (
    atom_id           BIGINT GENERATED BY DEFAULT AS IDENTITY PRIMARY KEY,
    structure_id      BIGINT NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    source_row        INTEGER,
    label             TEXT,
    element           TEXT,
    raw_type_symbol   TEXT,
    element_source    TEXT,
    x                 DOUBLE PRECISION,
    y                 DOUBLE PRECISION,
    z                 DOUBLE PRECISION,
    x_su              DOUBLE PRECISION,
    y_su              DOUBLE PRECISION,
    z_su              DOUBLE PRECISION,
    occupancy         DOUBLE PRECISION DEFAULT 1.0,
    occupancy_raw     TEXT,
    occupancy_state   TEXT NOT NULL DEFAULT 'legacy_unknown',
    u_iso             DOUBLE PRECISION,
    b_iso             DOUBLE PRECISION,
    adp_type          TEXT,
    disorder_assembly TEXT,
    disorder_group    TEXT,
    refinement_flags  TEXT,
    calc_flag         TEXT
);

-- Keep U/B data as supplied. atom_id is NULL when the label cannot be linked
-- unambiguously, but every row still belongs to its source structure.
CREATE TABLE atom_adp (
    adp_id          BIGINT GENERATED BY DEFAULT AS IDENTITY PRIMARY KEY,
    structure_id    BIGINT NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    atom_id         BIGINT REFERENCES atoms(atom_id) ON DELETE CASCADE,
    source_row      INTEGER,
    label           TEXT,
    representation  TEXT,
    u11             DOUBLE PRECISION,
    u22             DOUBLE PRECISION,
    u33             DOUBLE PRECISION,
    u12             DOUBLE PRECISION,
    u13             DOUBLE PRECISION,
    u23             DOUBLE PRECISION,
    b11             DOUBLE PRECISION,
    b22             DOUBLE PRECISION,
    b33             DOUBLE PRECISION,
    b12             DOUBLE PRECISION,
    b13             DOUBLE PRECISION,
    b23             DOUBLE PRECISION
);

CREATE TABLE import_issues (
    id          BIGSERIAL PRIMARY KEY,
    block_id    BIGINT NOT NULL REFERENCES cif_blocks(id) ON DELETE CASCADE,
    severity    TEXT NOT NULL CHECK (severity IN ('warning', 'error')),
    code        TEXT NOT NULL,
    source_row  INTEGER,
    tag         TEXT,
    message     TEXT NOT NULL
);

CREATE INDEX idx_cif_blocks_hash ON cif_blocks(sha256);
CREATE INDEX idx_import_issues_block ON import_issues(block_id);
CREATE INDEX idx_import_issues_code ON import_issues(code);

CREATE INDEX idx_structures_formula ON structures(formula);
CREATE INDEX idx_structures_space_group ON structures(space_group);
CREATE INDEX idx_structures_space_group_num ON structures(space_group_num);
CREATE INDEX idx_structures_centro ON structures(is_centrosymmetric);
CREATE INDEX idx_structures_volume ON structures(volume);
CREATE INDEX idx_structures_refcode ON structures(refcode);
CREATE UNIQUE INDEX idx_structures_source_block ON structures(source_block_id);
CREATE INDEX idx_structures_selection ON structures(is_centrosymmetric, r_factor_gt)
    WHERE extraction_status = 'ok';

CREATE INDEX idx_atoms_element ON atoms(element);
CREATE INDEX idx_atoms_structure ON atoms(structure_id);
CREATE UNIQUE INDEX idx_atoms_source_row ON atoms(structure_id, source_row);
CREATE INDEX idx_atoms_label ON atoms(structure_id, label);
CREATE UNIQUE INDEX idx_adp_source_row ON atom_adp(structure_id, source_row);

CREATE VIEW structures_summary AS
SELECT
    s.id,
    s.refcode,
    s.formula,
    s.space_group,
    s.space_group_num,
    s.point_group,
    s.crystal_system,
    s.is_centrosymmetric,
    s.a, s.b, s.c,
    s.alpha, s.beta, s.gamma,
    s.volume,
    s.reduced_volume,
    s.z,
    COUNT(a.atom_id) AS n_atoms
FROM structures s
LEFT JOIN atoms a ON a.structure_id = s.id
GROUP BY s.id;

INSERT INTO schema_migrations(version) VALUES (2);
COMMIT;
