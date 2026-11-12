"""Physical invariants, quality filters, sampling and CSV round-trip. No DB required."""
import hashlib
import importlib.util
import os
import json
from pathlib import Path
import tempfile
import unittest
import gemmi
import numpy as np
import pandas as pd
from unncdb_dataset import Calculation, Filters, InsufficientDataError, ReflectionGrid, Selection, StructureRejected, load_csv, save_csv
from unncdb_dataset.builder import build_records
from unncdb_dataset.diffraction import calculate_features
from unncdb_dataset.structure import prepare_structure


def cif(name="ABCDEF", group=1, sites="C1 C .13 .24 .35", extra="", atom_extra=(), a=10):
    operations = "\n".join("'" + op.triplet() + "'" for op in gemmi.find_spacegroup_by_number(group).operations())
    return f'''data_{name}
_cell_length_a {a}
_cell_length_b 11
_cell_length_c 12
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_refine_ls_R_factor_gt .03
_diffrn_radiation_probe x-ray
_symmetry_Int_Tables_number {group}
{extra}
loop_
_space_group_symop_operation_xyz
{operations}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{chr(10).join(atom_extra)}
{sites}
'''.encode()


def record(i, label, name=None):
    raw = cif(name or f"SAMPLE{i}", group=2 if label else 1, a=10+i/10)
    return ({"structure_id": i, "refcode": name or f"SAMPLE{i}",
             "is_centrosymmetric": bool(label), "space_group_num": 2 if label else 1,
             "cif_sha256": hashlib.sha256(raw).hexdigest()}, raw)


def small_dataset(records=None, n=4, **kwargs):
    records = records if records is not None else [record(i, i % 2) for i in range(1, 9)]
    return build_records(iter(records), selection=Selection(n=n, **kwargs), filters=Filters(),
                         grid=ReflectionGrid.cube(1), calculation=Calculation(anomalous=False))


class StructureTests(unittest.TestCase):
    def prepare(self, raw, **kwargs):
        return prepare_structure(raw, Filters(**kwargs), Calculation(anomalous=False))

    def assert_rejected(self, raw, code, **kwargs):
        with self.assertRaises(StructureRejected) as caught:
            self.prepare(raw, **kwargs)
        self.assertEqual(caught.exception.code, code)

    def test_periodic_duplicate_and_boundary(self):
        structure = self.prepare(cif(sites="C1 C 1.2 0 0\nC2 C .2 0 0\nO1 O -.1 0 0"))
        self.assertEqual(len(structure.coords), 2)
        self.assertEqual(structure.duplicates_removed, 1)
        self.assertTrue(((structure.coords >= 0) & (structure.coords < 1)).all())
        self.assert_rejected(cif(sites="C1 C .999999 0 0\nN1 N 0 0 0"), "mixed_site_collision")

    def test_collision_after_symmetry(self):
        self.assert_rejected(cif(group=2, sites="C1 C .2 .3 .4\nN1 N .8 .7 .6"), "mixed_site_collision")
        structure = self.prepare(cif(group=2, sites="C1 C 0 0 0"))
        self.assertEqual(len(structure.coords), 1)
        self.assertEqual(structure.duplicates_removed, 1)

    def test_quality_markers(self):
        self.assert_rejected(cif(sites="Ca1? Ca 0 0 0"), "ambiguous_atom_label")
        self.assert_rejected(cif(sites="C1 C ? 0 0"), "invalid_coordinates")
        self.assert_rejected(cif(sites="C1 C 0 0 0\nC1 C .2 0 0"), "duplicate_atom_label")
        self.assert_rejected(cif(extra="_chemical_name_common 'disordered compound'"), "disorder_description")
        self.assert_rejected(cif(sites="C1 C 0 0 0 1", atom_extra=("_atom_site_disorder_group",)), "disorder_marker")
        structure = self.prepare(cif(sites="C1 C 0 0 0 .", atom_extra=("_atom_site_disorder_group",)))
        self.assertEqual(len(structure.coords), 1)

    def test_occupancy(self):
        self.assert_rejected(cif(sites="C1 C 0 0 0 .5", atom_extra=("_atom_site_occupancy",)), "partial_occupancy")
        self.assert_rejected(cif(sites="C1 C 0 0 0 ?", atom_extra=("_atom_site_occupancy",)), "unknown_occupancy")
        structure = self.prepare(cif(sites="C1 C 0 0 0"))
        self.assertEqual(structure.occupancies[0], 1)
        self.assert_rejected(cif(sites="C1 C 0 0 0 1\nC2 C 1 0 0 .5", atom_extra=("_atom_site_occupancy",)),
                             "duplicate_occupancy_conflict", reject_partial_occupancy=False)

    def test_metadata_and_element_filters(self):
        self.assert_rejected(cif(), "required_elements", required_elements=("Ca",))
        self.assert_rejected(cif(), "excluded_elements", excluded_elements=("C",))
        with self.assertRaises(StructureRejected) as caught:
            prepare_structure(cif(), Filters(), Calculation(), {"is_centrosymmetric": True})
        self.assertEqual(caught.exception.code, "label_conflict")
        with self.assertRaises(StructureRejected) as caught:
            prepare_structure(cif(), Filters(), Calculation(), {"cif_sha256": "0"*64})
        self.assertEqual(caught.exception.code, "checksum_mismatch")

    def test_dense_duplicates_are_bounded(self):
        exact = "\n".join(f"C{i} C 0 0 0" for i in range(1500))
        structure = self.prepare(cif(sites=exact))
        self.assertEqual(len(structure.coords), 1)
        self.assertEqual(structure.duplicates_removed, 1499)
        close = "\n".join(f"C{i} C {i*1e-10} 0 0" for i in range(1500))
        self.assert_rejected(cif(sites=close), "dense_atom_collisions")

    def test_fingerprint_periodic_order_invariant(self):
        left = self.prepare(cif(sites="C1 C .2 .3 .4\nN1 N .6 .7 .8"))
        right = self.prepare(cif(sites="N2 N 1.6 -.3 .8\nC2 C 1.2 .3 .4"))
        self.assertEqual(left.fingerprint, right.fingerprint)


class PhysicsTests(unittest.TestCase):
    def test_nacl_extinctions_and_special_positions(self):
        structure = prepare_structure(cif(group=225, sites="Na1 Na 0 0 0\nCl1 Cl .5 0 0"), Filters(), Calculation(anomalous=False))
        self.assertEqual(len(structure.coords), 8)
        values = calculate_features(structure, ReflectionGrid([[1,0,0],[1,1,1],[2,0,0]]),
                                    Calculation(scattering="Z", anomalous=False, normalization="none"))
        np.testing.assert_allclose(values, [0,24**2,112**2], atol=1e-8)

    def test_friedel_pairs_real_and_anomalous(self):
        raw = cif(sites="C1 C .13 .24 .35\nBr1 Br .27 .38 .49")
        structure = prepare_structure(raw, Filters(), Calculation())
        grid = ReflectionGrid([[1,2,3],[-1,-2,-3]])
        real = calculate_features(structure, grid, Calculation(anomalous=False, normalization="none"))
        np.testing.assert_allclose(real[0], real[1], rtol=1e-12)
        anomalous = calculate_features(structure, grid, Calculation(normalization="none"))
        self.assertGreater(abs(anomalous[0]-anomalous[1]), 1e-4)

    def test_missing_reflection_is_not_zero(self):
        structure = prepare_structure(cif(), Filters(), Calculation())
        values = calculate_features(structure, ReflectionGrid([[1,0,0],[100,0,0]]), Calculation())
        self.assertEqual(values[0], 1)
        self.assertTrue(np.isnan(values[1]))

    def test_isotropic_damping(self):
        cfg = Calculation(anomalous=False, thermal="isotropic", normalization="none")
        raw = cif(sites="C1 C 0 0 0 .04", atom_extra=("_atom_site_U_iso_or_equiv",))
        structure = prepare_structure(raw, Filters(), cfg)
        grid = ReflectionGrid([[1,0,0],[2,0,0]])
        damped = calculate_features(structure, grid, cfg)
        ideal = calculate_features(structure, grid, Calculation(anomalous=False, normalization="none"))
        s2 = np.array([.05,.1])**2
        np.testing.assert_allclose(damped/ideal, np.exp(-16*np.pi**2*.04*s2), rtol=1e-12)

    @unittest.skipUnless(importlib.util.find_spec("pymatgen") and importlib.util.find_spec("tabulate"),
                         "Optional notebook comparison needs pymatgen and tabulate")
    def test_reference_notebook(self):
        notebook = Path(__file__).resolve().parents[2]/'src/str_amp.ipynb'
        ns = {}
        for cell in json.loads(notebook.read_text())["cells"][:7]:
            exec("".join(cell.get("source", [])), ns)
        raw = cif(sites="C1 C .13 .24 .35\nBr1 Br .27 .38 .49", group=2)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'reference.cif';path.write_bytes(raw)
            reference = ns['structural_amplitude_table'](path, hkl_max=2, theta_max=90,
                         wave_length=.71073, include_zero=False, unique_reflections=False)
        structure = prepare_structure(raw, Filters(), Calculation(anomalous=False))
        grid = ReflectionGrid.cube(2)
        values = calculate_features(structure, grid, Calculation(anomalous=False, normalization="none"))
        expected = {tuple(h): abs(f)**2 for h, f in zip(reference['hkl'], reference['F_hkl'])}
        np.testing.assert_allclose(values, [expected[tuple(h)] for h in grid.hkl], rtol=1e-12, atol=1e-8)


class DatasetTests(unittest.TestCase):
    def test_balance_features_and_repeatability(self):
        left, right = small_dataset(), small_dataset()
        self.assertIs(type(left), pd.DataFrame)
        pd.testing.assert_frame_equal(left, right)
        manifest = left.attrs["manifest"]
        features = left[manifest["feature_columns"]]
        self.assertEqual(left["is_centrosymmetric"].value_counts().to_dict(), {0: 2, 1: 2})
        self.assertEqual(features.shape, (4, 32))
        self.assertNotIn("space_group_num", features)
        self.assertNotIn("is_centrosymmetric", features)
        self.assertTrue(left[manifest["reflection_columns"]].notna().to_numpy().all())
        self.assertEqual(np.asarray(manifest["hkl"]).shape, (26, 3))

    def test_family_deduplication(self):
        records = [record(1,0,'ABCDEF01'),record(2,0,'ABCDEF02'),record(3,1,'GHIJKL'),record(4,0,'MNOPQR'),record(5,1,'STUVWX')]
        dataset = small_dataset(records)
        self.assertEqual(len(dataset["family_id"].unique()), 4)
        self.assertIn('duplicate_family', dataset.attrs["manifest"]["report"]['rejection_counts'])

    def test_partial_is_explicit(self):
        with self.assertRaises(InsufficientDataError) as caught:
            small_dataset([record(1,0)], n=4)
        self.assertFalse(caught.exception.report['complete'])
        self.assertEqual(len(caught.exception.partial_dataset), 1)
        dataset = small_dataset([record(1,0)], n=4, allow_partial=True)
        self.assertFalse(dataset.attrs["manifest"]["report"]['complete'])

    def test_csv_and_checksum(self):
        dataset = small_dataset()
        with tempfile.TemporaryDirectory() as d:
            path = save_csv(dataset, Path(d) / "dataset.csv")
            loaded = load_csv(path)
            manifest = dataset.attrs["manifest"]
            features = manifest["feature_columns"]
            self.assertIs(type(loaded), pd.DataFrame)
            self.assertEqual(loaded.attrs["manifest"]["hkl"], manifest["hkl"])
            np.testing.assert_allclose(dataset[features], loaded[features], rtol=0, atol=0)
            self.assertEqual(loaded["is_centrosymmetric"].tolist(), dataset["is_centrosymmetric"].tolist())
            self.assertEqual(loaded["family_id"].tolist(), dataset["family_id"].tolist())
            with self.assertRaises(FileExistsError):
                save_csv(dataset, path)
            path.write_text(path.read_text() + "\n")
            with self.assertRaises(ValueError):
                load_csv(path)

    def test_invalid_configuration(self):
        for value in (True,1.2,-1):
            with self.assertRaises(ValueError):ReflectionGrid.cube(value)
        with self.assertRaises(ValueError):ReflectionGrid([[0,0,0]])
        with self.assertRaises(ValueError):Selection(n=3)
        with self.assertRaises(ValueError):Calculation(scattering='Z')


@unittest.skipUnless(os.environ.get("UNNCDB_DATASET_TEST_DSN"), "Optional read-only database test")
class DatabaseTests(unittest.TestCase):
    def test_read_only_reproducible_sample(self):
        from unncdb_dataset import Database
        db = Database(os.environ["UNNCDB_DATASET_TEST_DSN"])
        with db._connection() as conn:
            self.assertEqual(conn.execute("SHOW transaction_read_only").fetchone()["transaction_read_only"], "on")
            self.assertEqual(conn.execute("SHOW transaction_isolation").fetchone()["transaction_isolation"], "repeatable read")
        audit = db.audit()
        self.assertGreater(int(audit.candidates.sum()), 0)
        first = db.build(n=4, grid=ReflectionGrid.cube(1), max_candidates=1000)
        second = db.build(n=4, grid=ReflectionGrid.cube(1), max_candidates=1000)
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(first["is_centrosymmetric"].value_counts().to_dict(), {0:2,1:2})


if __name__ == '__main__':
    unittest.main()
