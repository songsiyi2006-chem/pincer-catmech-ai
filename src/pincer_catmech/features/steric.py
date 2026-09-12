"""Molecular steric descriptors; Cartesian distances and radii are in angstrom.

Bondi ligand radii: A. Bondi, J. Phys. Chem. 68 (1964), 441-451,
doi:10.1021/j100785a001. The 1.17 scale, 3.5 A sphere, metal removal,
and optional hydrogens follow https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html.
This implementation integrates the union by Monte Carlo, not SambVca's mesh.
Sampling intervals quantify numerical error only, not chemical model uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from os import PathLike
from pathlib import Path
from statistics import NormalDist
from types import MappingProxyType
from typing import Mapping, Sequence

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from numpy.typing import ArrayLike, NDArray

# Deliberately limited to established ligand entries. No extended-table fallback.
BONDI_RADII: Mapping[str, float] = MappingProxyType({
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
})
Structure = Atoms | str | PathLike[str] | ArrayLike


@dataclass(frozen=True)
class BuriedVolumeResult:
    """Monte Carlo estimate; errors and interval endpoints are percentage points.

    The standard error is the plug-in binomial estimate and is zero at observed
    0/100%; the Wilson interval remains non-degenerate at these endpoints.
    Atom indices are zero-based and describe the atoms actually included.
    Effective radii include the scale, including for user-supplied overrides.
    """

    percent_buried_volume: float
    standard_error_percent: float
    confidence_interval_percent: tuple[float, float]
    confidence_level: float
    n_samples: int
    n_buried: int
    seed: int | None
    chunk_size: int
    metal_index: int
    atom_indices: tuple[int, ...]
    effective_radii: tuple[float, ...]
    sphere_radius: float
    radius_scale: float
    include_hydrogens: bool
    coordinate_units: str = "angstrom"
    sampling_method: str = "uniform sphere; NumPy PCG64"
    interval_method: str = "Wilson score"


def _positive(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive finite number")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive finite number") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number")
    return value


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _index(value: int, n_atoms: int) -> int:
    value = _integer(value, "atom index")
    if value >= n_atoms:
        raise ValueError(f"atom index {value} outside 0..{n_atoms - 1}")
    return value


def _geometry(
    structure: Structure, symbols: Sequence[str] | None, coordinate_units: str,
) -> tuple[NDArray[np.float64], tuple[str, ...]]:
    if coordinate_units != "angstrom":
        raise ValueError("coordinate_units must be 'angstrom'; convert coordinates explicitly")
    if isinstance(structure, Atoms):
        if symbols is not None:
            raise ValueError("symbols must be omitted for ASE Atoms")
        if np.any(structure.pbc):
            raise ValueError("periodic structures are unsupported; provide an isolated molecule")
        positions, symbols = structure.get_positions(), structure.get_chemical_symbols()
    elif isinstance(structure, (str, PathLike)):
        if symbols is not None:
            raise ValueError("symbols must be omitted for XYZ files")
        lines = Path(structure).read_text(encoding="utf-8-sig").splitlines()
        try:
            count = int(lines[0].strip())
            if count < 1 or len(lines) < count + 2 or any(x.strip() for x in lines[count + 2:]):
                raise ValueError("expected exactly one nonempty XYZ frame")
            rows = [line.split() for line in lines[2:count + 2]]
            if any(len(row) != 4 for row in rows):
                raise ValueError("XYZ rows require symbol and three Cartesian coordinates")
            symbols = [row[0] for row in rows]
            positions = [[float(x) for x in row[1:]] for row in rows]
        except (ValueError, IndexError) as exc:
            raise ValueError(f"malformed XYZ file: {exc}") from exc
    else:
        positions = structure
        if symbols is None:
            raise ValueError("raw Cartesian arrays require explicit symbols")
    try:
        positions = np.asarray(positions, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError("coordinates must be a numeric (N, 3) array") from exc
    if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) == 0:
        raise ValueError("coordinates must be a nonempty (N, 3) array")
    if not np.all(np.isfinite(positions)):
        raise ValueError("coordinates must be finite")
    if isinstance(symbols, str):
        raise ValueError("symbols must be a sequence of element symbols, not a formula")
    symbols = tuple(symbols)
    if len(symbols) != len(positions) or any(
        not isinstance(s, str) or s not in atomic_numbers or s == "X" for s in symbols
    ):
        raise ValueError("provide one valid element symbol per coordinate row")
    return positions, symbols


def bite_angle(
    structure: Structure, metal_index: int, donor1_index: int, donor2_index: int, *,
    symbols: Sequence[str] | None = None, coordinate_units: str = "angstrom",
) -> float:
    """Return the exact geometric L1-M-L2 angle in degrees; indices are zero-based.

    XYZ input must contain one ordinary frame (symbol x y z). Raw (N,3) arrays
    require symbols. Chemical bonding and atom identities are not inferred.
    """
    positions, _ = _geometry(structure, symbols, coordinate_units)
    metal, first, second = [_index(i, len(positions)) for i in (metal_index, donor1_index, donor2_index)]
    if len({metal, first, second}) != 3:
        raise ValueError("metal and both donor indices must be distinct")
    with np.errstate(over="ignore", invalid="ignore"):
        vectors = positions[[first, second]] - positions[metal]
    scales = np.max(np.abs(vectors), axis=1)
    if np.any(scales == 0) or not np.all(np.isfinite(scales)):
        raise ValueError("donor-metal vectors must be finite and nonzero")
    vectors = vectors / scales[:, None]
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    cosine = float(np.clip(vectors[0] @ vectors[1], -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def buried_volume(
    structure: Structure, metal_index: int, *, ligand_indices: Sequence[int],
    symbols: Sequence[str] | None = None, sphere_radius: float = 3.5,
    radius_scale: float = 1.17, include_hydrogens: bool = False,
    radii_overrides: Mapping[str, float] | None = None, n_samples: int = 100_000,
    seed: int | None = 0, chunk_size: int = 10_000, confidence_level: float = 0.95,
    coordinate_units: str = "angstrom",
) -> BuriedVolumeResult:
    """Estimate percent of the metal-centered sphere occupied by a selected ligand.

    Explicit ligand_indices prevent inclusion of solvent, counterions, or other
    ligands by accident. The central metal is always removed, even if selected.
    An empty selection is allowed and gives zero occupancy. Hydrogens are excluded
    by default. Unsupported included elements need an unscaled angstrom radius in
    radii_overrides; this includes any additional transition metal in the ligand.

    Uniform volume samples use z=2u-1, azimuth=2*pi*v, and r=R*w**(1/3).
    Sphere overlaps count once. Runtime is O(n_samples * n_ligand); temporary
    memory is O(chunk_size + n_ligand), not O(n_samples * n_ligand).
    Fixed seeds reproduce samples independently of chunk_size. No periodic
    imaging, conformational averaging, or wet-lab validation is implied.
    """
    positions, elements = _geometry(structure, symbols, coordinate_units)
    metal = _index(metal_index, len(positions))
    selected = tuple(_index(i, len(positions)) for i in ligand_indices)
    if len(selected) != len(set(selected)):
        raise ValueError("ligand_indices must not contain duplicates")
    if not isinstance(include_hydrogens, (bool, np.bool_)):
        raise ValueError("include_hydrogens must be boolean")
    selected = tuple(i for i in selected if i != metal and (include_hydrogens or elements[i] != "H"))
    sphere_radius = _positive(sphere_radius, "sphere_radius")
    radius_scale = _positive(radius_scale, "radius_scale")
    n_samples = _integer(n_samples, "n_samples", 1)
    chunk_size = _integer(chunk_size, "chunk_size", 1)
    seed = None if seed is None else _integer(seed, "seed")
    confidence_level = _positive(confidence_level, "confidence_level")
    if confidence_level >= 1:
        raise ValueError("confidence_level must be between zero and one")
    radii = dict(BONDI_RADII)
    for element, radius in (radii_overrides or {}).items():
        if element not in atomic_numbers or element == "X":
            raise ValueError(f"invalid radius override element: {element!r}")
        radii[element] = _positive(radius, f"radius for {element}")
    missing = sorted({elements[i] for i in selected} - radii.keys())
    if missing:
        raise ValueError(f"no supported Bondi radius for {', '.join(missing)}; provide radii_overrides")
    effective = np.array([radii[elements[i]] * radius_scale for i in selected])
    with np.errstate(over="ignore", invalid="ignore"):
        centers = positions[list(selected)] - positions[metal]
    if not np.all(np.isfinite(effective)) or not np.all(np.isfinite(centers)):
        raise ValueError("scaled radii and metal-relative coordinates must remain finite")
    rng = np.random.Generator(np.random.PCG64(seed))
    n_buried = 0
    for start in range(0, n_samples, chunk_size):
        count = min(chunk_size, n_samples - start)
        uniform = rng.random((count, 3))
        z = 2.0 * uniform[:, 0] - 1.0
        phi = 2.0 * np.pi * uniform[:, 1]
        radial = sphere_radius * np.cbrt(uniform[:, 2])
        planar = np.sqrt(np.maximum(0.0, 1.0 - z * z))
        points = radial[:, None] * np.column_stack((planar * np.cos(phi), planar * np.sin(phi), z))
        occupied = np.zeros(count, dtype=bool)
        for center, radius in zip(centers, effective):
            # One atom at a time bounds memory regardless of ligand atom count.
            with np.errstate(over="ignore"):
                occupied |= np.hypot.reduce(points - center, axis=1) <= radius
        n_buried += int(np.count_nonzero(occupied))
    proportion = n_buried / n_samples
    # Use lower-tail quantile to avoid rounding (1 + confidence_level)/2 to 1.
    quantile = -NormalDist().inv_cdf((1.0 - confidence_level) / 2.0)
    z2_n = quantile**2 / n_samples
    midpoint = (proportion + z2_n / 2.0) / (1.0 + z2_n)
    half_width = quantile * np.sqrt(proportion * (1.0 - proportion) / n_samples + z2_n / (4.0 * n_samples)) / (1.0 + z2_n)
    lower = 0.0 if n_buried == 0 else max(0.0, midpoint - half_width)
    upper = 1.0 if n_buried == n_samples else min(1.0, midpoint + half_width)
    return BuriedVolumeResult(
        percent_buried_volume=100.0 * proportion,
        standard_error_percent=float(100.0 * np.sqrt(proportion * (1.0 - proportion) / n_samples)),
        confidence_interval_percent=(100.0 * lower, 100.0 * upper),
        confidence_level=confidence_level, n_samples=n_samples, n_buried=n_buried,
        seed=seed, chunk_size=chunk_size, metal_index=metal, atom_indices=selected,
        effective_radii=tuple(float(x) for x in effective), sphere_radius=sphere_radius,
        radius_scale=radius_scale, include_hydrogens=bool(include_hydrogens),
    )
