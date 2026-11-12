"""Quality checks and periodic atom model preparation from authoritative CIF bytes."""
from dataclasses import dataclass
import hashlib
import json
import re
import gemmi
import numpy as np
from scipy.spatial import cKDTree
from .config import Calculation, Filters
from .errors import StructureRejected


@dataclass
class PreparedStructure:
    cell: gemmi.UnitCell
    symbols: np.ndarray
    coords: np.ndarray
    occupancies: np.ndarray
    u_iso: np.ndarray
    space_group_num: int
    centrosymmetric: bool
    duplicates_removed: int
    occupancy_assumptions: int
    fingerprint: str
    r_factor: float


def _reject(code, message):
    raise StructureRejected(code, message)


def _text(raw):
    return None if raw in (None, "", ".", "?") else gemmi.cif.as_string(raw)


def _number(raw, code="invalid_number"):
    if raw in (None, "", ".", "?"):
        _reject(code, "Missing numeric value")
    # Keep an exponent following a standard uncertainty, unlike truncation at '('.
    token = re.sub(r"\(\d+\)", "", raw)
    if not re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", token):
        _reject(code, f"Invalid numeric value: {raw!r}")
    value = float(token)
    if not np.isfinite(value):
        _reject(code, "Non-finite numeric value")
    return value


def _first(block, *tags):
    for tag in tags:
        raw = block.find_value(tag)
        if _text(raw) is not None:
            return raw
    return None


def _element(raw_type, label):
    raw = _text(raw_type) or label
    match = re.match(r"^([A-Z][a-z]?)", raw or "")
    if not match:
        _reject("unknown_element", f"Cannot determine element for {label!r}")
    symbol = match.group(1)
    if symbol in ("D", "T"):
        return "H"  # X-ray scattering; retain source isotope in raw CIF.
    element = gemmi.Element(symbol)
    if element.atomic_number == 0 or (raw_type and _text(raw_type) and not re.fullmatch(r"[A-Z][a-z]?[0-9+-]*", raw)):
        _reject("unknown_element", f"Unsupported atom type: {raw!r}")
    return element.name


def _operations(block, metadata):
    raw_ops = list(block.find_values("_space_group_symop_operation_xyz"))
    if not raw_ops:
        raw_ops = list(block.find_values("_symmetry_equiv_pos_as_xyz"))
    hall = _text(_first(block, "_space_group_name_Hall", "_symmetry_space_group_name_Hall"))
    try:
        if raw_ops:
            ops = [gemmi.Op(_text(value)).wrap() for value in raw_ops]
            unique = {op.triplet(): op for op in ops}
            ops = list(unique.values())
            group = gemmi.GroupOps(ops)
        elif hall:
            group = gemmi.symops_from_hall(hall)
            ops = list(group)
        else:
            _reject("missing_symmetry_operations", "Neither explicit operations nor Hall symbol is available")
        space_group = gemmi.find_spacegroup_by_ops(group)
        if space_group is None:
            _reject("unsupported_symmetry", "Operations do not form a supported complete group in this setting")
        if hall:
            hall_ops = {op.wrap().triplet() for op in gemmi.symops_from_hall(hall)}
            if hall_ops != {op.wrap().triplet() for op in ops}:
                _reject("symmetry_conflict", "Hall symbol and operations disagree")
        declared = _first(block, "_space_group_IT_number", "_symmetry_Int_Tables_number")
        if declared is not None and _number(declared) != space_group.number:
            _reject("symmetry_conflict", "Declared number and operations disagree")
        hm = _text(_first(block, "_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"))
        if hm:
            named = gemmi.find_spacegroup_by_name(hm)
            if named is not None and named.number != space_group.number:
                _reject("symmetry_conflict", "H-M symbol and operations disagree")
        centro = group.is_centrosymmetric()
        if metadata.get("space_group_num") is not None and int(metadata["space_group_num"]) != space_group.number:
            _reject("symmetry_conflict", "CIF and database group numbers disagree")
        if metadata.get("is_centrosymmetric") is not None and bool(metadata["is_centrosymmetric"]) != centro:
            _reject("label_conflict", "CIF and database target labels disagree")
        return ops, space_group.number, centro
    except (ValueError, RuntimeError, TypeError) as exc:
        if isinstance(exc, StructureRejected):
            raise
        _reject("invalid_symmetry", str(exc))


def prepare_structure(raw_cif: bytes | str, filters: Filters, calculation: Calculation,
                      metadata=None) -> PreparedStructure:
    metadata = metadata or {}
    raw_bytes = raw_cif.encode("utf-8") if isinstance(raw_cif, str) else bytes(raw_cif)
    if metadata.get("cif_sha256") and hashlib.sha256(raw_bytes).hexdigest() != metadata["cif_sha256"]:
        _reject("checksum_mismatch", "Stored CIF checksum does not match its bytes")
    try:
        document = gemmi.cif.read_string(raw_bytes.decode("utf-8"))
        if len(document) != 1:
            _reject("invalid_cif", "Expected exactly one CIF data block")
        block = document.sole_block()
    except (ValueError, RuntimeError, UnicodeError) as exc:
        if isinstance(exc, StructureRejected):
            raise
        _reject("invalid_cif", str(exc))
    r_factor = _number(block.find_value("_refine_ls_R_factor_gt"), "missing_r_factor")
    if not 0 <= r_factor <= filters.r_max:
        _reject("r_factor", f"R1_gt={r_factor} exceeds the permitted range")
    if metadata.get("r_factor_gt") is not None and not np.isclose(r_factor, metadata["r_factor_gt"], rtol=0, atol=1e-10):
        _reject("metadata_conflict", "CIF and database R factors disagree")
    if filters.radiation_probe:
        probe = _text(block.find_value("_diffrn_radiation_probe"))
        if not probe or probe.lower() != filters.radiation_probe.lower():
            _reject("radiation_probe", "Radiation probe does not match the selected experiment")
    if filters.reject_disorder:
        for tag in ("_chemical_name_systematic", "_chemical_name_common", "_refine_special_details",
                    "_publ_section_exptl_refinement"):
            value = _text(block.find_value(tag))
            if value and re.search(r"\bdisorder(?:ed|ing)?\b", value, flags=re.I):
                _reject("disorder_description", f"Disorder mentioned in {tag}")
    params = tuple(_number(block.find_value(tag), "invalid_cell") for tag in (
        "_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha",
        "_cell_angle_beta", "_cell_angle_gamma"))
    a, b, c, alpha, beta, gamma = params
    if min(a, b, c) <= 0 or not all(0 < v < 180 for v in (alpha, beta, gamma)):
        _reject("invalid_cell", "Cell lengths/angles are outside physical bounds")
    cosines = np.cos(np.radians([alpha, beta, gamma]))
    determinant = 1 + 2 * np.prod(cosines) - np.dot(cosines, cosines)
    if determinant <= 1e-12:
        _reject("invalid_cell", "Cell metric is singular or nonphysical")
    cell = gemmi.UnitCell(*params)
    if filters.max_volume is not None and cell.volume > filters.max_volume:
        _reject("cell_volume", "Cell exceeds max_volume")
    for name, value in zip(("a", "b", "c", "alpha", "beta", "gamma"), params):
        if metadata.get(name) is not None and not np.isclose(value, metadata[name], rtol=0, atol=1e-8):
            _reject("metadata_conflict", f"CIF and database {name} disagree")
    operations, group_num, centro = _operations(block, metadata)
    if filters.space_groups is not None and group_num not in filters.space_groups:
        _reject("space_group", "Space group is outside the requested set")
    tags = ("_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
            "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
            "_atom_site_disorder_group", "_atom_site_disorder_assembly", "_atom_site_calc_flag",
            "_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv")
    columns = {tag: list(block.find_values(tag)) for tag in tags}
    count = len(columns["_atom_site_label"])
    if not count:
        _reject("missing_atoms", "Atom-site labels are missing")
    if count > filters.max_asymmetric_atoms:
        _reject("atom_limit", "Independent atom table exceeds max_asymmetric_atoms")
    for tag, values in columns.items():
        if values and len(values) != count:
            _reject("atom_table_length", f"Column length mismatch: {tag}")
    for tag in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"):
        if len(columns[tag]) != count:
            _reject("missing_coordinates", f"Missing fractional coordinate column: {tag}")
    def value(tag, i):
        return columns[tag][i] if columns[tag] else None
    symbols, positions, occupancies, u_values = [], [], [], []
    labels = set()
    assumptions = 0
    for i in range(count):
        label = _text(value("_atom_site_label", i))
        if not label or "?" in label:
            _reject("ambiguous_atom_label", f"Missing/ambiguous label: {label!r}")
        if label in labels:
            _reject("duplicate_atom_label", f"Repeated atom label: {label}")
        labels.add(label)
        if filters.reject_disorder:
            for tag in ("_atom_site_disorder_group", "_atom_site_disorder_assembly"):
                marker = _text(value(tag, i))
                if marker not in (None, "0"):
                    _reject("disorder_marker", f"Disorder marker in {tag}: {marker}")
        if (_text(value("_atom_site_calc_flag", i)) or "").lower() == "dum":
            _reject("dummy_atom", "Dummy atom in source model")
        symbol = _element(value("_atom_site_type_symbol", i), label)
        position = [_number(value(f"_atom_site_fract_{axis}", i), "invalid_coordinates") for axis in "xyz"]
        raw_occupancy = value("_atom_site_occupancy", i)
        if raw_occupancy is None:
            occupancy = 1.0  # CIF dictionary default for absent column only.
        elif raw_occupancy in ("?", "."):
            if filters.require_known_occupancy:
                _reject("unknown_occupancy", "Explicit unknown/not-applicable occupancy")
            occupancy = 1.0
            assumptions += 1
        else:
            occupancy = _number(raw_occupancy, "invalid_occupancy")
        tol = filters.occupancy_tolerance
        if not -tol <= occupancy <= 1 + tol:
            _reject("invalid_occupancy", "Occupancy is outside [0,1] tolerance")
        if filters.reject_partial_occupancy and abs(occupancy - 1) > tol:
            _reject("partial_occupancy", "Atom site is not fully occupied")
        occupancy = float(np.clip(occupancy, 0, 1))
        u_iso = 0.0
        if calculation.thermal == "isotropic":
            raw_u, raw_b = value("_atom_site_U_iso_or_equiv", i), value("_atom_site_B_iso_or_equiv", i)
            if raw_u is not None:
                u_iso = _number(raw_u, "missing_thermal_parameters")
            elif raw_b is not None:
                u_iso = _number(raw_b, "missing_thermal_parameters") / (8 * np.pi**2)
            else:
                _reject("missing_thermal_parameters", "Isotropic mode requires U or B for every site")
            if u_iso < 0:
                _reject("invalid_thermal_parameters", "Negative isotropic displacement")
        symbols.append(symbol)
        positions.append(position)
        occupancies.append(occupancy)
        u_values.append(u_iso)
    present = set(symbols)
    if not set(filters.required_elements) <= present:
        _reject("required_elements", "Required elements are missing")
    if set(filters.excluded_elements) & present:
        _reject("excluded_elements", "Excluded elements are present")
    base_coords = np.asarray(positions, dtype=float)
    rotations = np.asarray([op.rot for op in operations], dtype=float) / gemmi.Op.DEN
    translations = np.asarray([op.tran for op in operations], dtype=float) / gemmi.Op.DEN
    # Bound temporary expansion too; maximum operation multiplicity is 192.
    if count * len(operations) > max(192 * filters.max_cell_atoms, 100000):
        _reject("atom_limit", "Symmetry expansion is too large")
    coords = np.mod(np.einsum("nj,sij->sni", base_coords, rotations) + translations[:, None, :], 1.0).reshape(-1, 3)
    coords[coords >= 1] = 0  # Floating-point modulo may round tiny negative values to 1.
    symbols = np.tile(np.asarray(symbols), len(operations))
    occupancies = np.tile(np.asarray(occupancies), len(operations))
    u_values = np.tile(np.asarray(u_values), len(operations))
    expanded_count = len(coords)
    # Collapse exact repetitions before enumerating close pairs. Special positions
    # and Mercury-expanded molecules can otherwise produce quadratic pair lists.
    _, first, inverse, multiplicities = np.unique(coords, axis=0, return_index=True,
                                                return_inverse=True, return_counts=True)
    if len(first) != len(coords):
        order = np.argsort(inverse, kind="stable")
        starts = np.r_[0, np.cumsum(multiplicities)[:-1]]
        _, symbol_ids = np.unique(symbols, return_inverse=True)
        if np.any(np.minimum.reduceat(symbol_ids[order], starts) != np.maximum.reduceat(symbol_ids[order], starts)):
            _reject("mixed_site_collision", "Different elements occupy the same periodic position")
        if np.any(np.maximum.reduceat(occupancies[order], starts) - np.minimum.reduceat(occupancies[order], starts) > filters.occupancy_tolerance):
            _reject("duplicate_occupancy_conflict", "Coincident sites have different occupancies")
        if calculation.thermal != "none" and np.any(np.maximum.reduceat(u_values[order], starts) - np.minimum.reduceat(u_values[order], starts) > 1e-8):
            _reject("duplicate_thermal_conflict", "Coincident sites have different displacement parameters")
        first = np.sort(first)
        symbols, coords, occupancies, u_values = symbols[first], coords[first], occupancies[first], u_values[first]
    tree = cKDTree(coords, boxsize=1.0)
    pair_count = (tree.count_neighbors(tree, filters.position_tolerance, p=np.inf) - len(coords)) // 2
    if pair_count > 1_000_000:
        _reject("dense_atom_collisions", "More than one million near-coincident atom pairs; unsafe/inconsistent model")
    pairs = tree.query_pairs(filters.position_tolerance, p=np.inf, output_type="ndarray")
    if len(pairs):
        if np.any(symbols[pairs[:, 0]] != symbols[pairs[:, 1]]):
            _reject("mixed_site_collision", "Different elements occupy the same periodic position")
        if np.any(np.abs(occupancies[pairs[:, 0]] - occupancies[pairs[:, 1]]) > filters.occupancy_tolerance):
            _reject("duplicate_occupancy_conflict", "Coincident sites have different occupancies")
        if calculation.thermal != "none" and np.any(np.abs(u_values[pairs[:, 0]] - u_values[pairs[:, 1]]) > 1e-8):
            _reject("duplicate_thermal_conflict", "Coincident sites have different displacement parameters")
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    keep = np.ones(len(coords), dtype=bool)
    for left, right in pairs:
        if keep[left]:
            keep[right] = False
    symbols, coords, occupancies, u_values = symbols[keep], coords[keep], occupancies[keep], u_values[keep]
    removed = expanded_count - len(coords)
    if len(coords) > filters.max_cell_atoms:
        _reject("atom_limit", "Expanded cell exceeds max_cell_atoms")
    # Exact-model deduplication, periodic coordinates/order insensitive. Not a general
    # crystallographic structure matcher: changed basis/origin and polymorphs need families.
    canonical = sorted((str(s), *np.mod(np.round(r, 8), 1.0).tolist(), round(float(o), 8),
                        round(float(u), 8)) for s, r, o, u in zip(symbols, coords, occupancies, u_values))
    fingerprint = hashlib.sha256(json.dumps([np.round(params, 6).tolist(), canonical],
                                           separators=(",", ":")).encode()).hexdigest()
    return PreparedStructure(cell, symbols, coords, occupancies, u_values, group_num, centro,
                             removed, assumptions, fingerprint, r_factor)
