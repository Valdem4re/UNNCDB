"""Validated, serializable settings. No model fitting or train/test splitting."""
from dataclasses import dataclass
import numbers
import numpy as np


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class Filters:
    r_max: float = 0.1
    reject_disorder: bool = True
    reject_partial_occupancy: bool = True
    require_known_occupancy: bool = True
    require_verified_symmetry: bool = True
    radiation_probe: str | None = "x-ray"
    space_groups: tuple[int, ...] | None = None
    required_elements: tuple[str, ...] = ()
    excluded_elements: tuple[str, ...] = ()
    position_tolerance: float = 1e-5  # Fractional, periodic L-infinity metric.
    occupancy_tolerance: float = 1e-6
    max_cell_atoms: int = 10000
    max_asymmetric_atoms: int = 5000
    max_volume: float | None = None

    def __post_init__(self):
        if not np.isfinite(self.r_max) or not 0 <= self.r_max <= 1:
            raise ValueError("r_max must be in [0, 1]")
        if not np.isfinite(self.position_tolerance) or not 0 < self.position_tolerance < 0.01:
            raise ValueError("position_tolerance must be in (0, 0.01)")
        if not np.isfinite(self.occupancy_tolerance) or not 0 <= self.occupancy_tolerance < 0.1:
            raise ValueError("occupancy_tolerance must be in [0, 0.1)")
        for name in ("max_cell_atoms", "max_asymmetric_atoms"):
            _positive_int(getattr(self, name), name)
        if self.max_volume is not None and (not np.isfinite(self.max_volume) or self.max_volume <= 0):
            raise ValueError("max_volume must be positive")
        if self.space_groups is not None:
            groups = tuple(self.space_groups)
            if not groups or any(isinstance(g, bool) or not isinstance(g, numbers.Integral)
                                 or not 1 <= g <= 230 for g in groups):
                raise ValueError("space_groups must contain integers in 1..230")
            object.__setattr__(self, "space_groups", groups)
        for name in ("required_elements", "excluded_elements"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        if set(self.required_elements) & set(self.excluded_elements):
            raise ValueError("required_elements and excluded_elements overlap")


@dataclass(frozen=True)
class Selection:
    n: int = 1000
    balance: bool = True
    seed: int = 42
    one_per_family: bool = True
    max_per_space_group: int | None = None
    max_candidates: int | None = None
    allow_partial: bool = False

    def __post_init__(self):
        _positive_int(self.n, "n")
        if self.balance and self.n % 2:
            raise ValueError("An exactly balanced sample requires even n")
        if isinstance(self.seed, bool) or not isinstance(self.seed, numbers.Integral) or not 0 <= self.seed < 2**63:
            raise ValueError("seed must be an integer in [0, 2**63)")
        for name in ("max_per_space_group", "max_candidates"):
            if getattr(self, name) is not None:
                _positive_int(getattr(self, name), name)


@dataclass(frozen=True)
class Calculation:
    wavelength: float = 0.71073  # Angstrom, synthetic experiment, not read from each record.
    scattering: str = "it92"
    anomalous: bool = True
    thermal: str = "none"  # Explicit ideal atoms, or complete isotropic U/B.
    quantity: str = "intensity"
    normalization: str = "mean_intensity"  # Not crystallographic E².
    theta_min: float = 0.0
    theta_max: float = 90.0
    batch_size: int = 128

    def __post_init__(self):
        if not np.isfinite(self.wavelength) or self.wavelength <= 0:
            raise ValueError("wavelength must be positive")
        if not np.isfinite([self.theta_min, self.theta_max]).all() or not 0 <= self.theta_min < self.theta_max <= 90:
            raise ValueError("Require 0 <= theta_min < theta_max <= 90 degrees")
        if self.scattering not in ("it92", "Z"):
            raise ValueError("scattering must be 'it92' or 'Z'")
        if self.scattering == "Z" and self.anomalous:
            raise ValueError("Z approximation requires anomalous=False")
        if self.thermal not in ("none", "isotropic"):
            raise ValueError("thermal must be 'none' or 'isotropic'")
        if self.quantity not in ("intensity", "amplitude"):
            raise ValueError("quantity must be 'intensity' or 'amplitude'")
        if self.normalization not in ("none", "mean_intensity"):
            raise ValueError("normalization must be 'none' or 'mean_intensity'")
        _positive_int(self.batch_size, "batch_size")


class ReflectionGrid:
    """One immutable hkl order shared by every structure; 000 is never a feature."""
    def __init__(self, hkl):
        values = np.asarray(hkl)
        if values.ndim != 2 or values.shape[1] != 3 or not len(values):
            raise ValueError("hkl must be a nonempty array of shape (N, 3)")
        if values.dtype.kind not in "iuf" or not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
            raise ValueError("hkl must contain finite integers")
        if np.any(values < -10000) or np.any(values > 10000) or len(values) > 100000:
            raise ValueError("Reflection grid is too large")
        values = values.astype(np.int32)
        if np.any(np.all(values == 0, axis=1)):
            raise ValueError("000 cannot be included in the feature grid")
        if len(np.unique(values, axis=0)) != len(values):
            raise ValueError("Repeated hkl in reflection grid")
        self._hkl = values.copy()
        self._hkl.flags.writeable = False

    @classmethod
    def cube(cls, hkl_max=4):
        try:
            values = (hkl_max,) * 3 if isinstance(hkl_max, numbers.Integral) else tuple(hkl_max)
        except TypeError as exc:
            raise ValueError("hkl_max must be an integer or three integers") from exc
        if len(values) != 3:
            raise ValueError("hkl_max must be an integer or three integers")
        for value in values:
            _positive_int(value, "hkl_max")
        if np.prod([2 * v + 1 for v in values], dtype=float) > 100001:
            raise ValueError("Reflection grid is too large")
        grid = np.stack(np.meshgrid(*(np.arange(-v, v + 1) for v in values), indexing="ij"), axis=-1).reshape(-1, 3)
        return cls(grid[np.any(grid != 0, axis=1)])

    @property
    def hkl(self):
        return self._hkl.copy()

    def __len__(self):
        return len(self._hkl)
