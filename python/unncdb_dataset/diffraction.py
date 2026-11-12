"""Fixed-grid kinematic X-ray scattering, based on the project's str_amp notebook."""
from functools import lru_cache
import gemmi
import numpy as np
from .config import Calculation, ReflectionGrid
from .errors import StructureRejected
from .structure import PreparedStructure


@lru_cache(maxsize=512)
def _coefficients(symbol, anomalous, wavelength):
    element = gemmi.Element(symbol)
    coeff = element.it92
    if coeff is None:
        raise StructureRejected("unsupported_scattering", f"No IT92 coefficients for {symbol}")
    fp, fpp = 0.0, 0.0
    if anomalous:
        z = element.atomic_number
        if z > 92:
            raise StructureRejected("unsupported_anomalous", "Cromer-Liberman tables support Z <= 92")
        if z >= 79 and not 0.16 <= wavelength <= 2.67:
            raise StructureRejected("unsupported_anomalous", "Wavelength is outside the supported heavy-element range")
        if z >= 3:
            fp, fpp = gemmi.cromer_liberman(z, gemmi.hc / wavelength)
    return np.asarray(coeff.a), np.asarray(coeff.b), coeff.c, complex(fp, fpp)


def calculate_features(structure: PreparedStructure, grid: ReflectionGrid, config: Calculation):
    hkl = grid.hkl
    # cell.orth.mat has lattice vectors as columns; reciprocal row vectors are h @ A^-1.
    stl = 0.5 * np.linalg.norm(hkl @ np.asarray(structure.cell.frac.mat), axis=1)
    low, high = np.sin(np.radians([config.theta_min, config.theta_max])) / config.wavelength
    observable = (stl >= low - 1e-12) & (stl <= high + 1e-12)
    selected = np.flatnonzero(observable)
    if not len(selected):
        raise StructureRejected("no_observable_reflections", "No grid reflection lies inside the experiment range")
    result = np.full(len(hkl), np.nan, dtype=np.float64)
    unique = np.unique(structure.symbols)
    # Remove cancellation noise relative to a forward-scattering amplitude bound.
    forward_bound = sum(gemmi.Element(str(s)).atomic_number * o for s, o in zip(structure.symbols, structure.occupancies))
    for start in range(0, len(selected), config.batch_size):
        ix = selected[start:start + config.batch_size]
        h = hkl[ix]
        s2 = stl[ix] ** 2
        total = np.zeros(len(ix), dtype=np.complex128)
        for symbol in unique:
            sites = structure.symbols == symbol
            if config.scattering == "Z":
                factors = np.full(len(ix), gemmi.Element(str(symbol)).atomic_number, dtype=complex)
            else:
                a, b, c, correction = _coefficients(str(symbol), config.anomalous, config.wavelength)
                factors = (a * np.exp(-s2[:, None] * b)).sum(axis=1) + c + correction
            phases = np.exp(2j * np.pi * (h @ structure.coords[sites].T))
            weights = structure.occupancies[sites]
            if config.thermal == "isotropic":
                phase_sum = (phases * np.exp(-8 * np.pi**2 * s2[:, None] * structure.u_iso[sites]) * weights).sum(axis=1)
            else:
                phase_sum = phases @ weights
            total += factors * phase_sum
        total[np.abs(total) <= max(forward_bound, 1.0) * 1e-12] = 0.0
        result[ix] = np.abs(total) ** 2
    if not np.isfinite(result[observable]).all():
        raise StructureRejected("nonfinite_features", "Non-finite calculated intensity")
    mean = result[observable].mean()
    if mean <= 0:
        raise StructureRejected("zero_signal", "All selected reflections have zero intensity")
    if config.normalization == "mean_intensity":
        result[observable] /= mean
    if config.quantity == "amplitude":
        result[observable] = np.sqrt(result[observable])
    return result
