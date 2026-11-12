"""Deterministic sampling after all physical/quality checks, with exact class quotas."""
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import importlib.metadata
import re
import time
import gemmi
import numpy as np
import pandas as pd
from .config import Calculation, Filters, ReflectionGrid, Selection
from .diffraction import calculate_features
from .errors import InsufficientDataError, StructureRejected
from .result import Dataset
from .structure import prepare_structure

CELL_COLUMNS = ["a", "b", "c", "alpha", "beta", "gamma"]
METADATA_COLUMNS = ["structure_id", "refcode", "family_id", "space_group_num", "is_centrosymmetric",
                    *CELL_COLUMNS, "r_factor_gt", "source_id", "source_block_id", "source_sha256",
                    "cif_sha256", "structure_hash", "atom_count", "duplicates_removed",
                    "occupancy_assumptions", "importer_version", "parser_version"]


@dataclass(frozen=True)
class BuildProgress:
    scanned: int
    accepted: int
    requested: int
    class_counts: dict
    elapsed_seconds: float


def family_id(refcode):
    """Conservative CSD refcode family, not a chemical similarity classifier."""
    name = str(refcode).upper()
    match = re.fullmatch(r"([A-Z]{6})(?:\d{2})?", name)
    return match.group(1) if match else name


def build_records(records, *, selection: Selection, filters: Filters, grid: ReflectionGrid,
                  calculation: Calculation, provenance=None, progress=None):
    start = time.monotonic()
    rows, vectors, rejections = [], [], []
    counts, group_counts = Counter(), Counter()
    families, fingerprints, source_hashes = set(), set(), set()
    scanned = 0
    def notify():
        if progress is not None:
            progress(BuildProgress(scanned, len(rows), selection.n,
                                   {0: counts[0], 1: counts[1]}, time.monotonic() - start))
    try:
        for metadata, raw in records:
            scanned += 1
            code = message = None
            label = metadata.get("is_centrosymmetric")
            family = family_id(metadata.get("refcode") or metadata.get("source_block_id"))
            if label is None:
                code, message = "missing_label", "No target label"
            elif selection.balance and counts[int(label)] >= selection.n // 2:
                code, message = "class_quota", "Class quota is already filled"
            elif selection.one_per_family and family in families:
                code, message = "duplicate_family", "Another accepted structure belongs to this refcode family"
            elif metadata.get("cif_sha256") and metadata["cif_sha256"] in source_hashes:
                code, message = "duplicate_cif", "CIF block is already represented"
            elif selection.max_per_space_group is not None and group_counts[metadata.get("space_group_num")] >= selection.max_per_space_group:
                code, message = "space_group_quota", "Space-group cap is filled"
            if code is None:
                try:
                    structure = prepare_structure(raw, filters, calculation, metadata)
                    if structure.fingerprint in fingerprints:
                        raise StructureRejected("duplicate_structure", "Same cell and periodic atomic model is already represented")
                    values = calculate_features(structure, grid, calculation)
                    label = int(structure.centrosymmetric)
                    row = {name: metadata.get(name) for name in METADATA_COLUMNS}
                    row.update(family_id=family, space_group_num=structure.space_group_num,
                               is_centrosymmetric=label, r_factor_gt=structure.r_factor,
                               structure_hash=structure.fingerprint, atom_count=len(structure.coords),
                               duplicates_removed=structure.duplicates_removed,
                               occupancy_assumptions=structure.occupancy_assumptions)
                    row.update(zip(CELL_COLUMNS, structure.cell.parameters))
                    rows.append(row)
                    vectors.append(values)
                    counts[label] += 1
                    group_counts[structure.space_group_num] += 1
                    families.add(family)
                    fingerprints.add(structure.fingerprint)
                    if metadata.get("cif_sha256"):
                        source_hashes.add(metadata["cif_sha256"])
                except StructureRejected as exc:
                    code, message = exc.code, str(exc)
            if code is not None:
                rejections.append({"structure_id": metadata.get("structure_id"),
                                   "refcode": metadata.get("refcode"), "code": code, "message": message})
            if scanned % 50 == 0:
                notify()
            if len(rows) == selection.n:
                break
    finally:
        close = getattr(records, "close", None)
        if close:
            close()
    notify()
    prefix = "I" if calculation.quantity == "intensity" else "F"
    reflection_columns = [f"{prefix}_{i:05d}" for i in range(len(grid))]
    matrix = np.asarray(vectors, dtype=np.float64).reshape(len(rows), len(grid))
    frame = pd.concat([pd.DataFrame(rows, columns=METADATA_COLUMNS),
                       pd.DataFrame(matrix, columns=reflection_columns)], axis=1)
    if rows:
        frame = frame.iloc[np.random.default_rng(selection.seed).permutation(len(frame))].reset_index(drop=True)
        frame["is_centrosymmetric"] = frame["is_centrosymmetric"].astype(np.uint8)
    complete = len(rows) == selection.n
    if selection.balance:
        complete = complete and counts[0] == counts[1]
    report = {"requested": selection.n, "accepted": len(rows), "scanned": scanned,
              "complete": complete, "class_counts": {str(k): counts[k] for k in (0, 1)},
              "rejection_counts": dict(Counter(r["code"] for r in rejections)),
              "elapsed_seconds": time.monotonic() - start}
    versions = {"numpy": np.__version__, "pandas": pd.__version__, "gemmi": gemmi.__version__,
                "unncdb_dataset": "0.1.0"}
    for name in ("scipy", "psycopg"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    manifest = {"format_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
                "filters": asdict(filters), "selection": asdict(selection),
                "calculation": asdict(calculation), "hkl": grid.hkl.tolist(),
                "reflection_columns": reflection_columns,
                "feature_columns": CELL_COLUMNS + reflection_columns,
                "target_column": "is_centrosymmetric", "target_meaning": {"0": "noncentrosymmetric", "1": "centrosymmetric"},
                "metadata_columns": METADATA_COLUMNS,
                "units": {"cell_lengths": "angstrom", "cell_angles": "degree", "wavelength": "angstrom",
                          "reflections": "dimensionless" if calculation.normalization != "none" else
                          "electrons_squared" if calculation.quantity == "intensity" else "electrons"},
                "missing_reflection": "NaN; not a zero intensity",
                "normalization_note": "mean_intensity divides I by the mean over observable fixed-grid reflections; not crystallographic E²",
                "symmetry_note": "Full grid, unmerged Friedel pairs, no true-group reflection reduction, 000 excluded",
                "family_note": "CSD base refcode; not a general chemical-family/structure matcher",
                "versions": versions, "database": provenance or {}, "report": report}
    result = Dataset(frame, manifest, pd.DataFrame(rejections, columns=["structure_id", "refcode", "code", "message"]))
    if not complete and not selection.allow_partial:
        exc = InsufficientDataError(f"Accepted {len(rows)}/{selection.n} structures after all checks; "
                                    f"class counts={dict(counts)}. Inspect exception.report or partial_dataset.", report)
        exc.partial_dataset = result
        raise exc
    return result
