"""Bounded CI-NEB/dimer searches with full-Hessian and reaction identity checks.

Every research energy/force comes from the injected ASE calculator. This module
does not infer ligand protonation sites, manufacture barriers, or accept a path
from its maximum-energy image alone. References: Henkelman et al., JCP 113 (2000)
9901, doi:10.1063/1.1329672; Henkelman/Jonsson, JCP 111 (1999) 7010,
doi:10.1063/1.480097; https://docs.ase-lib.org/ase/neb.html.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Callable, Mapping, Sequence

import numpy as np
from ase import Atoms, units
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate, NEB
from ase.optimize import FIRE

CalculatorFactory = Callable[[str], Calculator]
FREQUENCY_FACTOR = math.sqrt(units._e / (units._amu * 1e-20)) / (2 * math.pi * units._c * 100)
TRANSFER_KEYS = ("proton", "proton_donor", "proton_acceptor", "hydride", "hydride_donor", "hydride_acceptor")


class SearchDeadlineExceeded(TimeoutError):
    """A bounded search exhausted its wall-time allowance."""


@dataclass
class HessianResult:
    """Full Cartesian finite-difference Hessian and external-mode projection.

    Frequencies are signed cm^-1. Only rigid translation/rotation subspaces are
    removed; all internal atoms and modes remain. Cartesian eigenvectors have
    unit Euclidean norm, after undoing mass weighting.
    """

    hessian_ev_a2: np.ndarray
    frequencies_cm1: np.ndarray
    cartesian_modes: np.ndarray
    external_rank: int
    antisymmetry_relative: float
    displacement_angstrom: float

    def save(self, path: Path, atoms: Atoms) -> None:
        np.savez_compressed(path, hessian_ev_a2=self.hessian_ev_a2,
                            frequencies_cm1=self.frequencies_cm1, cartesian_modes=self.cartesian_modes,
                            external_rank=self.external_rank, antisymmetry_relative=self.antisymmetry_relative,
                            displacement_angstrom=self.displacement_angstrom,
                            positions=atoms.positions, atomic_numbers=atoms.numbers, masses=atoms.get_masses())


@dataclass
class PathSearchResult:
    status: str
    accepted: bool
    output_directory: str
    n_internal_images: int = 7
    neb_converged: bool = False
    provisional_neb_initial_guess: bool = False
    global_peak_image: int | None = None
    selected_initial_image: int | None = None
    restarted_ts_candidate: bool = False
    endpoint_converged: tuple[bool, bool] = (False, False)
    image_energies_ev: list[float] = field(default_factory=list)
    attempts: list[dict] = field(default_factory=list)
    failure_reasons: list[str] = field(default_factory=list)
    accepted_ts_energy_ev: float | None = None
    forward_electronic_barrier_ev: float | None = None
    reverse_electronic_barrier_ev: float | None = None
    elapsed_seconds: float = 0.0
    calculator_evaluations: int = 0
    evidence: str = "Injected-calculator electronic potential; no free-energy barrier inferred"


class _GuardedCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, backend: Calculator, guard: "_SearchState"):
        super().__init__()
        self.backend, self.guard = backend, guard

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        self.guard.check()
        super().calculate(atoms, properties, system_changes)
        self.backend.calculate(atoms, list(properties), system_changes)
        self.guard.evaluations += 1
        self.results = dict(self.backend.results)
        for name in ("energy", "forces"):
            if name in self.results and not np.all(np.isfinite(self.results[name])):
                raise ValueError(f"calculator returned nonfinite {name}")
        self.guard.check()


class _SearchState:
    def __init__(self, factory, directory, deadline_seconds, callback):
        self.factory, self.directory, self.callback = factory, directory, callback
        self.start = time.monotonic()
        self.deadline = math.inf if deadline_seconds is None else self.start + deadline_seconds
        self.evaluations = 0
        self.backends: list[Calculator] = []

    def check(self):
        if time.monotonic() >= self.deadline:
            raise SearchDeadlineExceeded("search wall-time allowance exhausted")

    def calculator(self, label):
        self.check()
        backend = self.factory(label)
        if not isinstance(backend, Calculator):
            raise TypeError("calculator_factory must return a fresh ASE Calculator")
        if any(backend is previous for previous in self.backends):
            raise ValueError("calculator_factory reused a calculator; each image needs independent state")
        self.backends.append(backend)
        return _GuardedCalculator(backend, self)

    def emit(self, stage, **values):
        self.check()
        record = {"stage": stage, "elapsed_seconds": time.monotonic() - self.start,
                  "calculator_evaluations": self.evaluations, **values}
        with (self.directory / "progress.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
        if self.callback is not None:
            self.callback(record)


def _molecule(atoms: Atoms) -> None:
    if not isinstance(atoms, Atoms) or len(atoms) < 2:
        raise ValueError("provide an ASE molecular Atoms object with at least two atoms")
    if np.any(atoms.pbc) or atoms.constraints:
        raise ValueError("full molecular Hessians require nonperiodic, unconstrained atoms")
    if not np.all(np.isfinite(atoms.positions)) or not np.all(np.isfinite(atoms.get_masses())):
        raise ValueError("coordinates and masses must be finite")
    if np.any(atoms.get_masses() <= 0):
        raise ValueError("atomic masses must be positive")


def full_hessian(atoms: Atoms, displacement: float = 0.01,
                 progress: Callable[[int, int], None] | None = None) -> HessianResult:
    """Central differences of true forces for every one of the 3N coordinates.

    Cost is 6N force evaluations. Energies use eV, distances angstrom, masses amu.
    The six rigid modes of nonlinear molecules (five for linear molecules) are
    projected in mass-weighted coordinates before diagonalization. No imaginary
    frequency is shifted, floored, removed, or converted to its absolute value.
    """
    _molecule(atoms)
    if not math.isfinite(displacement) or displacement <= 0:
        raise ValueError("displacement must be positive finite angstrom")
    if atoms.calc is None:
        raise ValueError("full_hessian requires a force-capable calculator")
    original = atoms.positions.copy()
    dimension = 3 * len(atoms)
    hessian = np.empty((dimension, dimension))
    try:
        for coordinate in range(dimension):
            moved = original.copy().reshape(-1)
            moved[coordinate] += displacement
            atoms.positions = moved.reshape(-1, 3)
            positive = atoms.get_forces().reshape(-1).copy()
            moved[coordinate] -= 2 * displacement
            atoms.positions = moved.reshape(-1, 3)
            negative = atoms.get_forces().reshape(-1).copy()
            hessian[:, coordinate] = -(positive - negative) / (2 * displacement)
            if progress is not None:
                progress(coordinate + 1, dimension)
    finally:
        atoms.positions = original
    return _project_cartesian_hessian(atoms, hessian, displacement)


def _project_cartesian_hessian(atoms, hessian, displacement):
    """Diagonalize the complete matrix; a provider cannot supply selected modes."""
    hessian = np.asarray(hessian, dtype=float)
    if hessian.shape != (3 * len(atoms), 3 * len(atoms)):
        raise ValueError("provider must supply a complete 3N by 3N Cartesian Hessian")
    if not np.all(np.isfinite(hessian)):
        raise ValueError("finite-difference Hessian contains nonfinite entries")
    asymmetry = float(np.linalg.norm(hessian - hessian.T) / max(np.linalg.norm(hessian), 1e-15))
    hessian = (hessian + hessian.T) / 2
    rootmass = np.sqrt(atoms.get_masses())
    repeated = np.repeat(rootmass, 3)
    weighted = hessian / np.outer(repeated, repeated)
    centered = atoms.positions - atoms.get_center_of_mass()
    external = []
    for axis in np.eye(3):
        external.append((rootmass[:, None] * np.broadcast_to(axis, (len(atoms), 3))).reshape(-1))
        external.append((rootmass[:, None] * np.cross(axis, centered)).reshape(-1))
    q, singular, _ = np.linalg.svd(np.column_stack(external), full_matrices=True)
    rank = int(np.sum(singular > singular[0] * 1e-9))
    if rank not in (5, 6):
        raise ValueError(f"invalid external molecular subspace rank {rank}; check geometry")
    internal = q[:, rank:]
    eigenvalues, eigenvectors = np.linalg.eigh(internal.T @ weighted @ internal)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * FREQUENCY_FACTOR
    modes = (internal @ eigenvectors) / repeated[:, None]
    modes /= np.linalg.norm(modes, axis=0)[None, :]
    return HessianResult(hessian, frequencies, modes.T.reshape(-1, len(atoms), 3),
                         rank, asymmetry, displacement)


def _transfer_spec(atoms: Atoms, indices: Mapping[str, int]) -> dict[str, int]:
    if set(indices) != set(TRANSFER_KEYS):
        raise ValueError(f"transfer_indices requires exactly {TRANSFER_KEYS}")
    result = {}
    for key, value in indices.items():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not 0 <= value < len(atoms):
            raise ValueError(f"invalid zero-based transfer index for {key}")
        result[key] = int(value)
    if result["proton"] == result["hydride"]:
        raise ValueError("proton and hydride must be different mapped atoms")
    for channel in ("proton", "hydride"):
        if atoms.numbers[result[channel]] != 1:
            raise ValueError(f"{channel} must reference a hydrogen atom")
        if len({result[channel], result[channel + "_donor"], result[channel + "_acceptor"]}) != 3:
            raise ValueError(f"{channel} donor, transferred atom and acceptor must be distinct")
    return result


def _progress_coordinate(atoms: Atoms, transfer: dict, channel: str) -> float:
    hydrogen = transfer[channel]
    return float(atoms.get_distance(hydrogen, transfer[channel + "_donor"])
                 - atoms.get_distance(hydrogen, transfer[channel + "_acceptor"]))


def transfer_mode_overlap(atoms: Atoms, mode: np.ndarray,
                          transfer_indices: Mapping[str, int]) -> dict[str, float | bool]:
    """Overlap with the Cartesian gradient of each donor-H minus acceptor-H distance.

    Both overlaps and their common direction are checked. The proton acceptor is
    explicitly supplied and may be N, O, or a dearomatizing C site.
    """
    transfer = _transfer_spec(atoms, transfer_indices)
    mode = np.asarray(mode, dtype=float)
    if mode.shape != (len(atoms), 3) or not np.all(np.isfinite(mode)) or np.linalg.norm(mode) == 0:
        raise ValueError("mode must be a nonzero finite (N,3) Cartesian displacement")
    mode = mode / np.linalg.norm(mode)
    overlaps = []
    gradients = []
    for channel in ("proton", "hydride"):
        gradient = np.zeros_like(mode)
        hydrogen = transfer[channel]
        for endpoint, sign in ((transfer[channel + "_donor"], 1), (transfer[channel + "_acceptor"], -1)):
            vector = atoms.positions[hydrogen] - atoms.positions[endpoint]
            length = np.linalg.norm(vector)
            if length < 1e-8:
                raise ValueError("overlapping atoms prevent transfer-coordinate definition")
            gradient[hydrogen] += sign * vector / length
            gradient[endpoint] -= sign * vector / length
        gradient /= np.linalg.norm(gradient)
        gradients.append(gradient)
        overlaps.append(float(np.sum(gradient * mode)))
    combined = gradients[0] + gradients[1]
    combined_overlap = abs(float(np.sum(combined * mode))) / np.linalg.norm(combined)
    return {"proton_overlap": abs(overlaps[0]), "hydride_overlap": abs(overlaps[1]),
            "combined_overlap": float(combined_overlap), "concerted_direction": overlaps[0] * overlaps[1] > 0}


def _aligned_rmsd(first: Atoms, second: Atoms) -> float:
    a, b = first.positions - first.positions.mean(axis=0), second.positions - second.positions.mean(axis=0)
    left, _, right = np.linalg.svd(a.T @ b)
    correction = np.diag([1.0, 1.0, np.linalg.det(left @ right)])
    return float(np.sqrt(np.mean(np.sum((a @ (left @ correction @ right) - b)**2, axis=1))))


def _optimize_minimum(atoms, directory, label, state, fmax, steps):
    atoms.calc = state.calculator(label)
    with FIRE(atoms, logfile=str(directory / f"{label}.log"), trajectory=str(directory / f"{label}.traj")) as optimizer:
        optimizer.attach(lambda: state.emit(label, step=optimizer.nsteps,
                                           energy_ev=float(atoms.get_potential_energy()),
                                           fmax_ev_a=float(np.linalg.norm(atoms.get_forces(), axis=1).max())))
        converged = bool(optimizer.run(fmax=fmax, steps=steps))
    write(directory / f"{label}.xyz", atoms)
    return converged


def refine_dimer(atoms: Atoms, initial_mode: np.ndarray, directory: Path,
                 fmax: float = 0.03, steps: int = 80) -> bool:
    """ASE minimum-mode following, using true forces and a supplied reaction tangent."""
    mode = np.asarray(initial_mode, dtype=float)
    if mode.shape != (len(atoms), 3) or not np.all(np.isfinite(mode)) or np.linalg.norm(mode) == 0:
        raise ValueError("initial dimer mode must be nonzero finite (N,3)")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with DimerControl(logfile=str(directory / "dimer_rotation.log"),
                      eigenmode_logfile=str(directory / "dimer_modes.log"),
                      f_rot_min=0.01, f_rot_max=0.05, max_num_rot=4,
                      dimer_separation=0.005, maximum_translation=0.05,
                      use_central_forces=False) as control:
        saddle = MinModeAtoms(atoms, control=control,
                              eigenmodes=[mode / np.linalg.norm(mode)], random_seed=0)
        with MinModeTranslate(saddle, logfile=str(directory / "dimer_translation.log"),
                              trajectory=str(directory / "dimer.traj")) as optimizer:
            converged = bool(optimizer.run(fmax=fmax, steps=steps))
    return converged and float(np.linalg.norm(atoms.get_forces(), axis=1).max()) <= fmax


def _descent_validation(ts, mode, endpoints, transfer, directory, state, fmax, steps,
                        displacement, rmsd_tolerance):
    branches = []
    for sign, label in ((-1, "minus"), (1, "plus")):
        branch = ts.copy()
        branch.positions += sign * displacement * mode / np.max(np.linalg.norm(mode, axis=1))
        converged = _optimize_minimum(branch, directory, f"{directory.name}_descent_{label}", state, fmax, steps)
        rmsds = [_aligned_rmsd(branch, endpoint) for endpoint in endpoints]
        progress = [_progress_coordinate(branch, transfer, channel) for channel in ("proton", "hydride")]
        nearest = int(np.argmin(rmsds))
        target = [_progress_coordinate(endpoints[nearest], transfer, channel) for channel in ("proton", "hydride")]
        branch_energy = float(branch.get_potential_energy())
        branches.append({"label": label, "converged": converged, "energy_ev": branch_energy,
                         "rmsd_to_endpoints_angstrom": rmsds, "nearest_endpoint": nearest,
                         "transfer_progress_angstrom": progress,
                         "same_transfer_basin": bool(all(a * b > 0 for a, b in zip(progress, target))),
                         "matches_endpoint": bool(rmsds[nearest] <= rmsd_tolerance),
                         "downhill": branch_energy < float(ts.get_potential_energy()) - 1e-5})
    accepted = all(branch["converged"] and branch["downhill"] and branch["same_transfer_basin"]
                   and branch["matches_endpoint"] for branch in branches)
    return accepted and {branch["nearest_endpoint"] for branch in branches} == {0, 1}, branches


def run_path(
    reactant: Atoms, product: Atoms, calculator_factory: CalculatorFactory,
    output_dir: str | Path, transfer_indices: Mapping[str, int], *,
    atom_mapping: Sequence[int] | None = None, n_internal_images: int = 7,
    neb_fmax: float = 0.05, neb_steps: int = 120, ts_fmax: float = 0.03,
    ts_steps: int = 80, endpoint_steps: int = 80, max_corrections: int = 2,
    hessian_step: float = 0.01, descent_steps: int = 80,
    descent_displacement: float = 0.15, endpoint_rmsd_tolerance: float = 0.5,
    imaginary_noise_cm1: float = 0.0, minimum_transfer_overlap: float = 0.15,
    minimum_combined_overlap: float = 0.30, deadline_seconds: float | None = None,
    progress_callback: Callable[[dict], None] | None = None,
    restart_images: Sequence[Atoms] | None = None,
    restart_source: str | Path | None = None,
    refine_unconverged_band: bool = False,
    ts_initial_image: int | None = None,
    restart_ts_candidate: Atoms | None = None,
    restart_ts_source: str | Path | None = None,
    restart_ts_mode: np.ndarray | None = None,
    restart_ts_mode_source: str | Path | None = None,
    hessian_provider: Callable[[Atoms, Path], HessianResult] | None = None,
) -> PathSearchResult:
    """Run seven-image CI-NEB and certify a concerted hydrogen-transfer saddle.

    ``atom_mapping[i]`` gives the product atom corresponding to reactant atom i;
    omission explicitly uses the supplied ordering. No element-based remapping
    is guessed. Endpoints must contain the same complete atoms/electronic state.
    Factory labels are unique; a deadline is checked before/after each backend
    evaluation (the backend must itself bound an individual SCF job).

    Strict default: every negative *internal* frequency counts as imaginary.
    A nonzero numerical-noise tolerance is available only as an explicit choice,
    recorded in the protocol. Acceptance additionally requires converged true
    forces, one imaginary in [-1800,-300] cm^-1, both transfer overlaps with the
    same sign, and converged downhill branches in opposite endpoint basins.
    Reported barriers, only on acceptance, are electronic eV, not qRRHO free
    energies. Numerical failures and exhausted searches are persisted as failures.
    A restart supplies all nine frames in reactant atom order, with first/last
    coordinates matching the supplied endpoints. Fresh calculators revalidate
    every frame; endpoint minimization exits immediately if fresh forces pass.
    Saved internal coordinates bypass interpolation and begin climbing NEB.
    An explicitly requested unconverged-band refinement uses its highest image
    only as a provisional dimer guess; NEB remains marked unconverged and every
    independent force/Hessian/reaction/descent saddle gate is still required.
    ts_initial_image explicitly selects an internal image (1..7) as the dimer
    guess; the global peak remains separately recorded. This can target a local
    chemical step in a band containing more than one structural rearrangement.
    A saved dimer candidate can be supplied with a restart band and its source
    file. Its actual coordinates are retained, while dimer optimizer history
    and eigenmode are reinitialized from the band tangent and fresh forces.
    restart_ts_mode optionally supplies an explicit Cartesian (N,3) direction
    for the first dimer attempt at a saved candidate. Its source file and both
    original/normalized mode arrays are archived. This changes initialization
    only; it never substitutes for the final full Hessian or reaction gates.
    An optional Hessian provider may use a native complete numerical Hessian.
    It must bound its own subprocess, preserve the candidate, return all atoms
    and modes in HessianResult, and archive native provenance in attempt_dir.
    The supplied matrix is independently projected again before TS validation.
    """
    _molecule(reactant)
    _molecule(product)
    if hessian_provider is not None and not callable(hessian_provider):
        raise ValueError("hessian_provider must be callable")
    if not isinstance(refine_unconverged_band, bool):
        raise ValueError("refine_unconverged_band must be boolean")
    if ts_initial_image is not None and (isinstance(ts_initial_image, bool) or
            not isinstance(ts_initial_image, int) or not 1 <= ts_initial_image <= 7):
        raise ValueError("ts_initial_image must be an internal image index from 1 through 7")
    for name, value in (("neb_steps", neb_steps), ("ts_steps", ts_steps), ("endpoint_steps", endpoint_steps),
                        ("descent_steps", descent_steps), ("max_corrections", max_corrections)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if n_internal_images != 7:
        raise ValueError("campaign protocol requires seven internal images (nine including endpoints)")
    for name, value in (("neb_fmax", neb_fmax), ("ts_fmax", ts_fmax), ("hessian_step", hessian_step),
                        ("descent_displacement", descent_displacement), ("endpoint_rmsd_tolerance", endpoint_rmsd_tolerance)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if imaginary_noise_cm1 < 0 or not math.isfinite(imaginary_noise_cm1):
        raise ValueError("imaginary_noise_cm1 must be nonnegative and finite")
    if not all(0 <= value <= 1 for value in (minimum_transfer_overlap, minimum_combined_overlap)):
        raise ValueError("overlap thresholds must lie in [0,1]")
    if deadline_seconds is not None and (not math.isfinite(deadline_seconds) or deadline_seconds <= 0):
        raise ValueError("deadline_seconds must be positive and finite")
    endpoints = [reactant.copy(), product.copy()]
    mapping = None
    if atom_mapping is not None:
        mapping = list(atom_mapping)
        if any(isinstance(i, bool) or not isinstance(i, (int, np.integer)) for i in mapping) or sorted(mapping) != list(range(len(product))):
            raise ValueError("atom_mapping must be a complete permutation of product indices")
        endpoints[1] = endpoints[1][mapping]
    if not np.array_equal(endpoints[0].numbers, endpoints[1].numbers) or not np.array_equal(endpoints[0].get_masses(), endpoints[1].get_masses()):
        raise ValueError("mapped endpoints must have identical atomic identities and masses")
    for key in ("charge", "unpaired", "uhf"):
        if endpoints[0].info.get(key) != endpoints[1].info.get(key):
            raise ValueError(f"endpoint electronic-state metadata differs: {key}")
    transfer = _transfer_spec(endpoints[0], transfer_indices)
    restart = None
    restart_provenance = None
    if restart_images is not None:
        restart = list(restart_images)
        if len(restart) != 9:
            raise ValueError("restart_images must contain exactly nine full molecular frames")
        for index, frame in enumerate(restart):
            _molecule(frame)
            if not np.array_equal(frame.numbers, endpoints[0].numbers) or not np.array_equal(frame.get_masses(), endpoints[0].get_masses()):
                raise ValueError(f"restart frame {index} has different atom identities, order, or masses")
            for key in ("charge", "unpaired", "uhf"):
                if frame.info.get(key) != endpoints[0].info.get(key):
                    raise ValueError(f"restart frame {index} electronic-state metadata differs: {key}")
        for frame, endpoint in ((restart[0], endpoints[0]), (restart[-1], endpoints[1])):
            if not np.allclose(frame.positions, endpoint.positions, rtol=0, atol=1e-8):
                raise ValueError("restart endpoint coordinates must match the supplied mapped endpoints")
        restart = [frame.copy() for frame in restart]
        restart_provenance = {"frames": 9, "atom_order": "reactant order", "interpolation": "none",
                              "initial_climbing_image": True, "all_energies_and_forces_recomputed": True}
        if restart_source is not None:
            source_path = Path(restart_source).resolve()
            restart_provenance.update(source=str(source_path), source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest())
    elif restart_source is not None:
        raise ValueError("restart_source requires restart_images")
    candidate_provenance = None
    if restart_ts_candidate is not None:
        if restart is None or restart_ts_source is None:
            raise ValueError("restart_ts_candidate requires restart_images and restart_ts_source")
        _molecule(restart_ts_candidate)
        if not np.array_equal(restart_ts_candidate.numbers, endpoints[0].numbers) or not np.array_equal(restart_ts_candidate.get_masses(), endpoints[0].get_masses()):
            raise ValueError("restart TS candidate atom identities or masses differ")
        for key in ("charge", "unpaired", "uhf"):
            if restart_ts_candidate.info.get(key) != endpoints[0].info.get(key):
                raise ValueError(f"restart TS candidate electronic-state metadata differs: {key}")
        source_path = Path(restart_ts_source).resolve()
        candidate_provenance = {"source": str(source_path), "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                                "optimizer_history": "reinitialized", "initial_eigenmode": "saved band tangent", "fresh_force_validation": True}
    elif restart_ts_source is not None:
        raise ValueError("restart_ts_source requires restart_ts_candidate")
    initial_mode = None
    mode_provenance = None
    if restart_ts_mode is not None:
        if restart_ts_candidate is None or restart_ts_mode_source is None:
            raise ValueError("restart_ts_mode requires restart_ts_candidate and restart_ts_mode_source")
        input_mode = np.asarray(restart_ts_mode, dtype=float).copy()
        if input_mode.shape != (len(endpoints[0]), 3) or not np.all(np.isfinite(input_mode)):
            raise ValueError("restart_ts_mode must be a finite Cartesian (N,3) array")
        scale = float(np.max(np.abs(input_mode)))
        if scale == 0:
            raise ValueError("restart_ts_mode must be nonzero")
        initial_mode = input_mode / scale
        initial_mode /= np.linalg.norm(initial_mode)
        mode_source = Path(restart_ts_mode_source).resolve()
        mode_provenance = {"source": str(mode_source),
                           "source_sha256": hashlib.sha256(mode_source.read_bytes()).hexdigest(),
                           "coordinates": "Cartesian atom order of restart candidate",
                           "normalization": "unit Euclidean norm over all 3N coordinates",
                           "usage": "first dimer attempt only; independent final saddle validation unchanged"}
        candidate_provenance["initial_eigenmode"] = "explicit saved Cartesian mode"
    elif restart_ts_mode_source is not None:
        raise ValueError("restart_ts_mode_source requires restart_ts_mode")
    directory = Path(output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / "result.json").exists() or (directory / "protocol.json").exists():
        raise FileExistsError("output directory contains a previous search; choose a new directory")
    provider_name = (f"{getattr(hessian_provider, '__module__', type(hessian_provider).__module__)}."
                     f"{getattr(hessian_provider, '__qualname__', type(hessian_provider).__qualname__)}") if hessian_provider else None
    protocol = {"method": "ASE CI-NEB plus dimer plus full force-difference Hessian", "n_internal_images": 7,
                "hessian_provider": provider_name or "ASE central differences of all Cartesian true forces",
                "transfer_indices": transfer, "atom_mapping": None if mapping is None else list(map(int, mapping)),
                "neb_fmax_ev_a": neb_fmax, "neb_steps": neb_steps, "ts_fmax_ev_a": ts_fmax,
                "ts_steps": ts_steps, "endpoint_steps": endpoint_steps, "max_corrections": max_corrections,
                "hessian_step_angstrom": hessian_step, "imaginary_noise_cm1": imaginary_noise_cm1,
                "imaginary_window_cm1": [-1800, -300], "descent_steps": descent_steps,
                "descent_displacement_angstrom": descent_displacement,
                "endpoint_rmsd_tolerance_angstrom": endpoint_rmsd_tolerance,
                "minimum_transfer_overlap": minimum_transfer_overlap, "minimum_combined_overlap": minimum_combined_overlap,
                "deadline_seconds": deadline_seconds, "free_energy_barrier_computed": False,
                "refine_unconverged_band": refine_unconverged_band,
                "requested_ts_initial_image": ts_initial_image,
                "restart_ts_candidate": candidate_provenance,
                "restart_ts_mode": mode_provenance,
                "restart": restart_provenance}
    if restart is not None:
        snapshot = directory / "restart_input.traj"
        write(snapshot, restart)
        restart_provenance.update(snapshot=snapshot.name, snapshot_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest())
    if restart_ts_candidate is not None:
        snapshot = directory / "restart_ts_input.traj"
        write(snapshot, restart_ts_candidate.copy())
        candidate_provenance.update(snapshot=snapshot.name, snapshot_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest())
    if initial_mode is not None:
        snapshot = directory / "restart_ts_mode.npz"
        np.savez_compressed(snapshot, input_mode=input_mode, normalized_mode=initial_mode,
                            positions=restart_ts_candidate.positions, atomic_numbers=restart_ts_candidate.numbers,
                            masses=restart_ts_candidate.get_masses())
        mode_provenance.update(snapshot=snapshot.name, snapshot_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest())
    (directory / "protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")
    result = PathSearchResult("started", False, str(directory))
    state = _SearchState(calculator_factory, directory, deadline_seconds, progress_callback)
    try:
        for index, endpoint in enumerate(endpoints):
            write(directory / f"endpoint_{index}_input.xyz", endpoint)
        endpoint_status = []
        for index, endpoint in enumerate(endpoints):
            endpoint_status.append(_optimize_minimum(endpoint, directory, f"endpoint_{index}", state, neb_fmax, endpoint_steps))
        result.endpoint_converged = tuple(endpoint_status)
        if not all(endpoint_status):
            raise ValueError("endpoint optimization did not converge within the step allowance")
        for channel in ("proton", "hydride"):
            if not (_progress_coordinate(endpoints[0], transfer, channel) < -0.10
                    and _progress_coordinate(endpoints[1], transfer, channel) > 0.10):
                raise ValueError(f"optimized endpoints do not bracket the annotated {channel} transfer")
        images = [endpoints[0]] + ([endpoints[0].copy() for _ in range(7)] if restart is None else restart[1:-1]) + [endpoints[1]]
        for index, image in enumerate(images[1:-1], 1):
            image.calc = state.calculator(f"neb_image_{index:02d}")
        neb = NEB(images, climb=restart is not None, parallel=False, method="improvedtangent",
                  remove_rotation_and_translation=True, allow_shared_calculator=False)
        if restart is None:
            neb.interpolate(method="idpp", apply_constraint=False)
        warmup_steps = min(20, neb_steps // 4)
        with FIRE(neb, logfile=str(directory / "neb.log"), trajectory=str(directory / "neb.traj")) as optimizer:
            optimizer.attach(lambda: state.emit("neb", step=optimizer.nsteps,
                                               image_energies_ev=[float(image.get_potential_energy()) for image in images]))
            if restart is None:
                optimizer.run(fmax=neb_fmax, steps=warmup_steps)
            neb.climb = True
            result.neb_converged = bool(optimizer.run(fmax=neb_fmax, steps=neb_steps - optimizer.nsteps))
        result.image_energies_ev = [float(image.get_potential_energy()) for image in images]
        write(directory / "neb_final.extxyz", images)
        if not result.neb_converged and not refine_unconverged_band:
            raise ValueError("CI-NEB did not converge within the step allowance")
        result.provisional_neb_initial_guess = not result.neb_converged
        if result.provisional_neb_initial_guess:
            state.emit("unconverged_band_provisional_guess", neb_converged=False,
                       independent_full_saddle_validation_required=True)
        result.global_peak_image = int(np.argmax(result.image_energies_ev[1:-1])) + 1
        peak = result.global_peak_image if ts_initial_image is None else ts_initial_image
        result.selected_initial_image = peak if restart_ts_candidate is None else None
        result.restarted_ts_candidate = restart_ts_candidate is not None
        state.emit("dimer_guess_selected", global_peak_image=result.global_peak_image,
                   global_peak_energy_ev=result.image_energies_ev[result.global_peak_image],
                   selected_image=result.selected_initial_image,
                   selected_energy_ev=None if restart_ts_candidate is not None else result.image_energies_ev[peak],
                   initial_tangent_image=peak,
                   selection=("saved dimer candidate" if restart_ts_candidate is not None else
                              "global maximum" if ts_initial_image is None else "explicit local image"))
        tangent = images[peak + 1].positions - images[peak - 1].positions
        if np.linalg.norm(tangent) < 1e-8:
            raise ValueError("vanishing reaction tangent at candidate saddle")
        candidate = images[peak].copy() if restart_ts_candidate is None else restart_ts_candidate.copy()
        correction_mode = None
        for attempt_index in range(max_corrections + 1):
            state.check()
            attempt_dir = directory / f"attempt_{attempt_index:02d}"
            attempt_dir.mkdir()
            attempt = {"attempt": attempt_index, "accepted": False, "status": "started",
                       "initial_guess_from_unconverged_band": result.provisional_neb_initial_guess,
                       "selected_initial_image": result.selected_initial_image, "global_peak_image": result.global_peak_image,
                       "restarted_ts_candidate": result.restarted_ts_candidate}
            result.attempts.append(attempt)
            candidate.calc = state.calculator(f"ts_attempt_{attempt_index:02d}")
            if correction_mode is not None:
                # Alternate signs without selecting a fabricated outcome; revalidate every new saddle.
                candidate.positions += (-1 if attempt_index % 2 else 1) * 0.12 * correction_mode
                attempt["orthogonal_correction_angstrom"] = 0.12
            state.emit("dimer_start", attempt=attempt_index)
            dimer_mode = initial_mode if attempt_index == 0 and initial_mode is not None else tangent
            attempt["initial_eigenmode"] = "explicit saved Cartesian mode" if dimer_mode is initial_mode else "band tangent"
            refined = refine_dimer(candidate, dimer_mode, attempt_dir, ts_fmax, ts_steps)
            attempt["dimer_converged"] = refined
            attempt["true_fmax_ev_a"] = float(np.linalg.norm(candidate.get_forces(), axis=1).max())
            attempt["energy_ev"] = float(candidate.get_potential_energy())
            state.emit("dimer_complete", attempt=attempt_index, converged=refined,
                       energy_ev=attempt["energy_ev"], fmax_ev_a=attempt["true_fmax_ev_a"])
            write(attempt_dir / "candidate.xyz", candidate)
            if not refined or attempt["true_fmax_ev_a"] > ts_fmax:
                attempt["status"] = "dimer_not_converged"
                continue
            if hessian_provider is None:
                hessian = full_hessian(candidate, hessian_step,
                    lambda done, total: state.emit("hessian", attempt=attempt_index, completed_columns=done, total_columns=total))
            else:
                state.emit("hessian_provider_start", attempt=attempt_index, provider=provider_name)
                supplied = candidate.copy()
                supplied.calc = candidate.calc
                hessian = hessian_provider(supplied, attempt_dir)
                state.check()
                if not (np.array_equal(supplied.positions, candidate.positions) and
                        np.array_equal(supplied.numbers, candidate.numbers) and
                        np.array_equal(supplied.get_masses(), candidate.get_masses()) and
                        all(supplied.info.get(key) == candidate.info.get(key) for key in ("charge", "unpaired", "uhf"))):
                    raise ValueError("Hessian provider changed the candidate geometry or atom identities")
                if not isinstance(hessian, HessianResult):
                    raise TypeError("Hessian provider must return HessianResult")
                if not (math.isfinite(hessian.displacement_angstrom) and hessian.displacement_angstrom > 0 and
                        math.isfinite(hessian.antisymmetry_relative) and hessian.antisymmetry_relative >= 0):
                    raise ValueError("invalid provider displacement or Hessian antisymmetry")
                verified = _project_cartesian_hessian(candidate, hessian.hessian_ev_a2, hessian.displacement_angstrom)
                if (hessian.external_rank != verified.external_rank or
                        np.shape(hessian.frequencies_cm1) != verified.frequencies_cm1.shape or
                        not np.allclose(hessian.frequencies_cm1, verified.frequencies_cm1, atol=1e-4, rtol=1e-7)):
                    raise ValueError("provider spectrum disagrees with the complete Cartesian Hessian")
                verified.antisymmetry_relative = max(verified.antisymmetry_relative, hessian.antisymmetry_relative)
                hessian = verified
                state.emit("hessian_provider_complete", attempt=attempt_index, provider=provider_name)
            hessian.save(attempt_dir / "full_hessian.npz", candidate)
            attempt.update(hessian_provider=provider_name or "ASE central differences",
                           hessian_displacement_angstrom=hessian.displacement_angstrom,
                           hessian_file="full_hessian.npz",
                           hessian_sha256=hashlib.sha256((attempt_dir / "full_hessian.npz").read_bytes()).hexdigest())
            frequencies = hessian.frequencies_cm1
            imaginary = np.flatnonzero(frequencies < -imaginary_noise_cm1)
            attempt.update(frequencies_cm1=frequencies.tolist(), imaginary_count=int(len(imaginary)),
                           raw_negative_count=int(np.sum(frequencies < 0)), external_rank=hessian.external_rank,
                           hessian_antisymmetry_relative=hessian.antisymmetry_relative)
            reaction_index = int(np.argmax([transfer_mode_overlap(candidate, mode, transfer)["combined_overlap"]
                                            for mode in hessian.cartesian_modes]))
            correction_indices = [i for i in imaginary if i != reaction_index]
            correction_mode = None
            if correction_indices:
                correction_mode = hessian.cartesian_modes[correction_indices[0]].copy()
                unit_tangent = tangent / np.linalg.norm(tangent)
                correction_mode -= np.sum(correction_mode * unit_tangent) * unit_tangent
                length = np.linalg.norm(correction_mode)
                correction_mode = correction_mode / length if length > 1e-8 else None
            valid_hessian = (len(imaginary) == 1 and -1800 <= frequencies[imaginary[0]] <= -300
                             and np.all(frequencies != 0) and hessian.antisymmetry_relative <= 0.05)
            if not valid_hessian:
                attempt["status"] = "hessian_rejected"
                continue
            mode = hessian.cartesian_modes[imaginary[0]]
            overlap = transfer_mode_overlap(candidate, mode, transfer)
            attempt["transfer_mode"] = overlap
            if not (overlap["proton_overlap"] >= minimum_transfer_overlap and overlap["hydride_overlap"] >= minimum_transfer_overlap
                    and overlap["combined_overlap"] >= minimum_combined_overlap and overlap["concerted_direction"]):
                attempt["status"] = "wrong_imaginary_mode"
                continue
            branches_ok, branches = _descent_validation(candidate, mode, endpoints, transfer, attempt_dir, state,
                                                         neb_fmax, descent_steps, descent_displacement, endpoint_rmsd_tolerance)
            attempt["descent_branches"] = branches
            if not branches_ok:
                attempt["status"] = "endpoint_descent_rejected"
                continue
            endpoint_energies = [float(endpoint.get_potential_energy()) for endpoint in endpoints]
            if attempt["energy_ev"] <= max(endpoint_energies):
                attempt["status"] = "saddle_not_above_endpoints"
                continue
            attempt.update(status="accepted", accepted=True)
            result.accepted = True
            result.accepted_ts_energy_ev = attempt["energy_ev"]
            result.forward_electronic_barrier_ev = attempt["energy_ev"] - endpoint_energies[0]
            result.reverse_electronic_barrier_ev = attempt["energy_ev"] - endpoint_energies[1]
            write(directory / "accepted_ts.xyz", candidate)
            break
        result.status = "accepted" if result.accepted else "rejected_after_bounded_corrections"
        if not result.accepted:
            result.failure_reasons.extend(attempt["status"] for attempt in result.attempts)
    except (SearchDeadlineExceeded, TimeoutError) as exc:
        result.status = "deadline_exhausted"
        result.failure_reasons.append(f"{type(exc).__name__}: {exc}")
        if result.attempts and result.attempts[-1]["status"] == "started":
            result.attempts[-1].update(status="deadline_exhausted", error=str(exc))
    except Exception as exc:
        result.status = "failed"
        result.failure_reasons.append(f"{type(exc).__name__}: {exc}")
        if result.attempts and result.attempts[-1]["status"] == "started":
            result.attempts[-1].update(status="failed", error=f"{type(exc).__name__}: {exc}")
    finally:
        result.elapsed_seconds = time.monotonic() - state.start
        result.calculator_evaluations = state.evaluations
        (directory / "result.json").write_text(json.dumps(asdict(result), indent=2, allow_nan=False), encoding="utf-8")
    return result
