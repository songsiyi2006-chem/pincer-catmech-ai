"""Native GFN2-xTB Hessians and explicit molecular RRHO/qRRHO thermochemistry.

The Hessian is computed by xTB --hess, not by free_energy.py. The unprojected
DFTB-style hessian.out is requested explicitly, read in Hartree/bohr^2, then
mass weighted and projected onto the full molecular internal subspace.
ALPB 'gsolv' contributes excess solvation only; the 1 atm -> 1 M ideal standard
state correction is subsequently added once by calculate_thermochemistry.
Sources: https://xtb-docs.readthedocs.io/en/latest/hessian.html and
https://xtb-docs.readthedocs.io/en/latest/gbsa.html#reference-states.
"""

from __future__ import annotations

from dataclasses import asdict
import gzip
import hashlib
import inspect
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Mapping, Sequence

import numpy as np
from ase import Atoms, units
from ase.io import read, write
from ase.thermochemistry import IdealGasThermo

from pincer_catmech.kinetics.free_energy import HARTREE_TO_KCAL, calculate_thermochemistry
from pincer_catmech.kinetics.neb_ts_search import FREQUENCY_FACTOR, HessianResult, _molecule
from pincer_catmech.quantum.xtb_backend import (
    NUMBER, XTBCalculator, atomic_json, calculation_environment, executable_path,
    optimize_geometry, parse_energy, utc_now, validate_electrons,
)

EV_TO_KCAL = HARTREE_TO_KCAL / units.Hartree
DEFAULT_TEMPERATURES = (298.15, 340.0, 350.0, 360.0, 370.0, 380.0, 383.15,
                        390.0, 400.0, 410.0, 420.0, 430.0, 440.0)


def read_native_hessian(path: str | Path, n_atoms: int) -> np.ndarray:
    """Read xTB's unprojected DFTB-style full Hessian; return eV/angstrom^2."""
    if isinstance(n_atoms, bool) or not isinstance(n_atoms, int) or n_atoms < 1:
        raise ValueError("n_atoms must be a positive integer")
    tokens = Path(path).read_text(encoding="utf-8").split()
    if len(tokens) != (3 * n_atoms)**2 or any(re.fullmatch(NUMBER, token) is None for token in tokens):
        raise ValueError("native unprojected Hessian must contain exactly (3N)^2 finite numeric entries")
    matrix = np.array([float(token.replace("D", "E").replace("d", "e")) for token in tokens])
    if not np.all(np.isfinite(matrix)):
        raise ValueError("native Hessian contains nonfinite entries")
    return matrix.reshape(3 * n_atoms, 3 * n_atoms) * units.Hartree / units.Bohr**2


def _rigid_basis(atoms):
    """Use one rank criterion for projection and RRHO mode cardinality."""
    rootmass = np.sqrt(atoms.get_masses())
    repeated = np.repeat(rootmass, 3)
    centered = atoms.positions - atoms.get_center_of_mass()
    external = []
    for axis in np.eye(3):
        external.extend([(rootmass[:, None] * np.broadcast_to(axis, (len(atoms), 3))).ravel(),
                         (rootmass[:, None] * np.cross(axis, centered)).ravel()])
    basis, singular, _ = np.linalg.svd(np.column_stack(external), full_matrices=True)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-9))
    if rank not in (5, 6):
        raise ValueError(f"invalid rigid molecular subspace rank: {rank}")
    return basis, rank, repeated


def project_hessian(atoms: Atoms, hessian_ev_a2: np.ndarray) -> HessianResult:
    """Mass-weight full Cartesian Hessian and remove only rigid external modes."""
    _molecule(atoms)
    matrix = np.asarray(hessian_ev_a2, dtype=float)
    dimension = 3 * len(atoms)
    if matrix.shape != (dimension, dimension) or not np.all(np.isfinite(matrix)):
        raise ValueError("full finite (3N,3N) Cartesian Hessian required")
    asymmetry = float(np.linalg.norm(matrix - matrix.T) / max(np.linalg.norm(matrix), 1e-15))
    symmetric = (matrix + matrix.T) / 2
    basis, rank, repeated = _rigid_basis(atoms)
    internal = basis[:, rank:]
    weighted = symmetric / np.outer(repeated, repeated)
    values, vectors = np.linalg.eigh(internal.T @ weighted @ internal)
    frequencies = np.sign(values) * np.sqrt(np.abs(values)) * FREQUENCY_FACTOR
    modes = (internal @ vectors) / repeated[:, None]
    modes /= np.linalg.norm(modes, axis=0)[None, :]
    # The native finite-difference step is recorded separately in native input.
    return HessianResult(symmetric, frequencies, modes.T.reshape(-1, len(atoms), 3), rank, asymmetry, 0.005 * units.Bohr)


def thermochemistry_from_modes(
    atoms: Atoms, energy_ev: float, frequencies_cm1: Sequence[float], *,
    transition_state: bool = False, unpaired: int = 0, symmetry_number: int = 1,
    temperatures: Sequence[float] = DEFAULT_TEMPERATURES, pressure_pa: float = 101325.0,
    concentration_mol_l: float = 1.0, cutoff_cm1: float = 100.0,
    source: str = "native GFN2-xTB Hessian", source_sha256: str | None = None,
    solvent: str | None = "toluene", solvation_state: str = "gsolv",
    ts_imaginary_window_cm1: tuple[float, float] | None = (-1800.0, -300.0),
) -> list:
    """Full molecular H/S plus entropy-only Grimme correction at each actual T.

    The electronic spin degeneracy is the explicit nominal 2S+1=unpaired+1,
    not a computed ground-spin assignment. Symmetry number one is an explicit
    default assumption. Transition-state mode lists must include their one
    negative mode, which is excluded from the partition function exactly once.
    None explicitly disables only the TS frequency window; independent chemical
    mode and descent validation remains the caller's responsibility.
    """
    _molecule(atoms)
    if not isinstance(transition_state, bool):
        raise ValueError("transition_state must be boolean")
    if isinstance(unpaired, bool) or not isinstance(unpaired, int) or unpaired < 0:
        raise ValueError("unpaired must be a nonnegative integer")
    if isinstance(symmetry_number, bool) or not isinstance(symmetry_number, int) or symmetry_number < 1:
        raise ValueError("symmetry_number must be a positive integer")
    if not math.isfinite(energy_ev):
        raise ValueError("electronic energy must be finite eV")
    if solvation_state != "gsolv":
        raise ValueError("use ALPB gsolv; other states already include a concentration shift")
    frequencies = np.asarray(frequencies_cm1, dtype=float)
    if frequencies.ndim != 1 or not np.all(np.isfinite(frequencies)):
        raise ValueError("frequencies must be a finite one-dimensional full internal spectrum")
    _, external_rank, _ = _rigid_basis(atoms)
    linear = external_rank == 5
    expected = 3 * len(atoms) - (5 if linear else 6)
    if len(frequencies) != expected:
        raise ValueError(f"full internal spectrum requires {expected} modes; received {len(frequencies)}")
    imaginary = frequencies[frequencies < 0]
    if len(imaginary) != int(transition_state) or np.any(frequencies == 0):
        raise ValueError("minimum needs all positive modes; TS needs exactly one negative and no zero modes")
    if ts_imaginary_window_cm1 is not None:
        window = np.asarray(ts_imaginary_window_cm1, dtype=float)
        if window.shape != (2,) or not np.all(np.isfinite(window)) or not window[0] < window[1] < 0:
            raise ValueError("TS imaginary window must contain two ordered, finite negative cm^-1 limits")
        ts_imaginary_window_cm1 = tuple(map(float, window))
        if transition_state and not window[0] <= imaginary[0] <= window[1]:
            raise ValueError(f"TS imaginary frequency must lie in {tuple(window)} cm^-1")
    positive = frequencies[frequencies > 0]
    kwargs = dict(vib_energies=(positive * units.invcm).tolist(),
                  geometry="linear" if linear else "nonlinear", potentialenergy=float(energy_ev),
                  atoms=atoms, symmetrynumber=symmetry_number, spin=unpaired / 2, ignore_imag_modes=False)
    # New ASE supports explicit all-mode selection; older ASE leaves a shorter
    # supplied TS list intact. Cardinality has already been checked above.
    if "vib_selection" in inspect.signature(IdealGasThermo).parameters:
        kwargs["vib_selection"] = "all"
    harmonic = IdealGasThermo(**kwargs)
    notes = ["Electronic spin occupation is assumed, not an independently determined ground state",
             f"Rotational symmetry number {symmetry_number} supplied explicitly",
             "ALPB excess solvation contribution reused across temperature; no temperature-dependent solvent model calculated"]
    if solvent is None:
        notes[-1] = "Gas-phase electronic potential; no solvent stabilization included"
    if transition_state:
        notes.append(f"TS imaginary frequency window: {ts_imaginary_window_cm1}; chemical mode and descent validation must be supplied independently")
    output = []
    seen = set()
    for temperature in temperatures:
        temperature = float(temperature)
        if not math.isfinite(temperature) or temperature <= 0 or temperature in seen:
            raise ValueError("temperature grid must contain distinct positive finite kelvin values")
        seen.add(temperature)
        enthalpy = harmonic.get_enthalpy(temperature, verbose=False) * EV_TO_KCAL
        entropy = harmonic.get_entropy(temperature, pressure=pressure_pa, verbose=False) * EV_TO_KCAL * 1000
        output.append(calculate_thermochemistry(
            frequencies.tolist(), E_elec=energy_ev * EV_TO_KCAL,
            ZPVE=harmonic.get_ZPE_correction() * EV_TO_KCAL, H_298=enthalpy,
            S_harmonic=entropy, temperature=temperature, pressure_pa=pressure_pa,
            concentration_mol_l=concentration_mol_l, cutoff_cm1=cutoff_cm1,
            transition_state=transition_state, source=source, source_sha256=source_sha256,
            program="GFN2-xTB+ALPB" if solvent else "GFN2-xTB",
            entropy_provenance="ASE full molecular RRHO translation/rotation/vibration/electronic entropy",
            warnings=notes,
        ))
    if not output:
        raise ValueError("temperature grid must not be empty")
    return output


def _remaining(deadline_epoch, limit):
    if deadline_epoch is None:
        return limit
    remaining = deadline_epoch - time.time()
    if remaining <= 0:
        raise TimeoutError("campaign deadline reached")
    return min(limit, remaining)


def _archive(path: Path, output: Path) -> dict:
    data = path.read_bytes()
    with gzip.open(output, "wb") as stream:
        stream.write(data)
    return {"file": output.name, "source_sha256": hashlib.sha256(data).hexdigest(),
            "source_bytes": len(data), "gzip_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}


def native_hessian(atoms, scratch_dir, *, charge=0, unpaired=0, threads=1,
                   executable=None, deadline_epoch=None, timeout_seconds=1800,
                   solvent="toluene", solvation_state="gsolv"):
    """Execute full unbiased xTB --hess with native unprojected matrix output."""
    _molecule(atoms)
    validate_electrons(atoms, charge, unpaired)
    if solvation_state != "gsolv":
        raise ValueError("only gsolv avoids a duplicate concentration correction")
    folder = Path(scratch_dir)
    folder.mkdir(parents=True, exist_ok=False)
    executable = executable_path(executable)
    write(folder / "input.xyz", atoms, format="xyz")
    (folder / "hessian.inp").write_text(
        "$hess\n  sccacc=0.05\n  step=0.005\n  scale=1.0\n"
        "$write\n  hessian.out=true\n$end\n", encoding="ascii")
    command = [executable, "input.xyz", "--gfn", "2", "--chrg", str(charge), "--uhf", str(unpaired),
               "--etemp", "300", "--iterations", "500", "--acc", "0.05", "--hess", "--input", "hessian.inp"]
    if solvent:
        command += ["--alpb", solvent, solvation_state]
    started = time.monotonic()
    atomic_json(folder / "command.json", {"argv": command, "started_utc": utc_now(),
                                         "threads": threads, "solvent": solvent, "solvation_state": solvation_state})
    try:
        completed = subprocess.run(command, cwd=folder,
            env=calculation_environment(executable, threads), capture_output=True,
            text=True, encoding="utf-8", errors="replace",
            timeout=_remaining(deadline_epoch, timeout_seconds))
        output = completed.stdout + "\n" + completed.stderr
    except subprocess.TimeoutExpired as exc:
        parts = [value.encode("utf-8") if isinstance(value, str) else value or b""
                 for value in (exc.stdout, exc.stderr)]
        (folder / "xtb.out").write_bytes(b"\n".join(parts))
        raise TimeoutError("native Hessian exceeded its wall-time allowance") from exc
    (folder / "xtb.out").write_text(output, encoding="utf-8")
    if completed.returncode or "normal termination of xtb" not in output.lower():
        raise RuntimeError("native xTB Hessian did not terminate normally")
    if re.search(r"biased numerical hessian|O1NumHess", output, re.I):
        raise ValueError("biased or O(1) approximate Hessian is not the requested full unbiased Hessian")
    frozen = re.search(rf"frozen atoms in %\s*:\s*({NUMBER})\s+(\d+)", output, re.I)
    scale = re.search(rf"Hessian scale factor\s*:\s*({NUMBER})", output, re.I)
    if frozen is None or int(frozen[2]) != 0 or float(frozen[1]) != 0:
        raise ValueError("native output must confirm zero frozen atoms for the full molecular Hessian")
    if scale is None or not math.isclose(float(scale[1]), 1.0, abs_tol=1e-8):
        raise ValueError("native output must confirm an unscaled Hessian")
    matrix = read_native_hessian(folder / "hessian.out", len(atoms))
    return project_hessian(atoms, matrix), parse_energy(output), time.monotonic() - started


def _coordination_check(atoms, original, specification):
    if specification is None:
        return {"checked": False, "retained": None}
    metal = int(specification["metal_index"])
    donors = [int(i) for i in specification["donor_indices"]]
    distances = [float(atoms.get_distance(metal, donor)) for donor in donors]
    baseline = [float(original.get_distance(metal, donor)) for donor in donors]
    retained = all(1.1 < current < max(3.0, 1.25 * prior) for current, prior in zip(distances, baseline))
    return {"checked": True, "retained": retained, "metal_donor_distances_A": distances,
            "criterion": "Geometric retention screen: 1.1 < d < max(3.0,1.25*d_input) angstrom; not a bonding proof"}


def characterize_stationary_point(
    atoms: Atoms, output_dir: str | Path, scratch_dir: str | Path, *,
    charge: int = 0, unpaired: int = 0, transition_state: bool = False,
    max_repairs: int = 2, temperatures: Sequence[float] = DEFAULT_TEMPERATURES,
    threads: int = 1, executable: str | None = None, deadline_epoch: float | None = None,
    solvent: str | None = "toluene", solvation_state: str = "gsolv",
    force_threshold_ev_a: float = 0.03, symmetry_number: int = 1,
    coordination_indices: Mapping | None = None, source_metadata: Mapping | None = None,
    ts_imaginary_window_cm1: tuple[float, float] | None = (-1800.0, -300.0),
    identity_metadata: Mapping | None = None,
) -> dict:
    """Certify a full molecular minimum/saddle before emitting thermochemistry.

    All negative internal modes remain in the record. Small imaginary modes are
    numerical_uncertainty, never silently ignored. At most two native reoptimiza-
    tions can follow mode displacements for minima; TS geometries are not repaired
    by minimization. First-order-saddle status is spectral only: callers must also
    require the NEB transfer-mode and endpoint-descent acceptance record.
    """
    _molecule(atoms)
    validate_electrons(atoms, charge, unpaired)
    if isinstance(max_repairs, bool) or not isinstance(max_repairs, int) or not 0 <= max_repairs <= 2:
        raise ValueError("max_repairs must be an integer from zero through two")
    if solvation_state != "gsolv":
        raise ValueError("gsolv is required before the one-time ideal concentration correction")
    if not math.isfinite(force_threshold_ev_a) or force_threshold_ev_a <= 0:
        raise ValueError("force threshold must be positive finite eV/angstrom")
    if ts_imaginary_window_cm1 is not None:
        window = np.asarray(ts_imaginary_window_cm1, dtype=float)
        if window.shape != (2,) or not np.all(np.isfinite(window)) or not window[0] < window[1] < 0:
            raise ValueError("TS imaginary window must contain two ordered, finite negative cm^-1 limits")
        ts_imaginary_window_cm1 = tuple(map(float, window))
    directory, scratch = Path(output_dir), Path(scratch_dir)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / "result.json").exists():
        raise FileExistsError("preserve previous stationary-point records; choose a fresh output directory")
    scratch.mkdir(parents=True, exist_ok=True)
    original, current = atoms.copy(), atoms.copy()
    write(directory / "input.xyz", original)
    started = time.monotonic()
    result = {"status": "started", "accepted": False, "method": "GFN2-xTB",
              "solvent": solvent, "solvation_state": solvation_state, "electronic_temperature_K": 300.0,
              "charge": charge, "unpaired_electrons": unpaired,
              "spin_provenance": "Nominal occupation assumption; no ground-spin determination",
              "symmetry_number": symmetry_number, "symmetry_provenance": "Explicit assumed symmetry number",
              "transition_state": transition_state, "force_threshold_eV_A": force_threshold_ev_a,
              "ts_imaginary_window_cm1": ts_imaginary_window_cm1,
              "temperatures_K": list(temperatures), "attempts": [], "thermochemistry": [],
              "source_metadata": dict(source_metadata or {}), "started_utc": utc_now(),
              "chemical_identity_metadata": dict(identity_metadata) if identity_metadata is not None else None,
              "standard_state": "ALPB gsolv excess contribution plus ideal 1 atm to 1 M correction once",
              "frequency_unit": "cm^-1", "energy_unit": "kcal/mol", "entropy_unit": "cal/(mol K)",
              "hessian_method": "Native full unbiased numerical xTB --hess; unprojected hessian.out mass weighted and projected",
              "input_geometry_sha256": hashlib.sha256((directory / "input.xyz").read_bytes()).hexdigest()}
    try:
        for attempt_index in range(1 + (0 if transition_state else max_repairs)):
            _remaining(deadline_epoch, 1)
            attempt = {"attempt": attempt_index, "status": "started"}
            result["attempts"].append(attempt)
            atomic_json(directory / "result.json", result)
            if identity_metadata is not None:
                from pincer_catmech.quantum.identity import catalyst_identity_screen
                identity = catalyst_identity_screen(current, identity_metadata)
                attempt["chemical_identity"] = identity
                if not identity["retained"]:
                    attempt["status"] = result["status"] = "chemical_identity_not_retained"
                    break
            current.calc = XTBCalculator(scratch / f"gradient_{attempt_index:02d}", charge=charge, unpaired=unpaired,
                threads=threads, executable=executable, deadline=deadline_epoch,
                solvent=solvent, solvation_state=solvation_state,
                trace_path=directory / f"gradient_{attempt_index:02d}.jsonl.gz")
            force = float(np.linalg.norm(current.get_forces(), axis=1).max())
            energy_gradient = float(current.get_potential_energy())
            attempt.update(max_force_eV_A=force, gradient_energy_eV=energy_gradient)
            attempt["phase"] = "native_hessian"
            atomic_json(directory / "result.json", result)
            hessian_folder = scratch / f"native_hessian_{attempt_index:02d}"
            spectrum, energy, wall = native_hessian(current, hessian_folder, charge=charge, unpaired=unpaired,
                threads=threads, executable=executable, deadline_epoch=deadline_epoch,
                solvent=solvent, solvation_state=solvation_state)
            spectrum_path = directory / f"hessian_{attempt_index:02d}.npz"
            spectrum.save(spectrum_path, current)
            attempt["spectrum_npz_sha256"] = hashlib.sha256(spectrum_path.read_bytes()).hexdigest()
            for source, stem in (("hessian.out", "unprojected_hessian"), ("xtb.out", "native_output"), ("hessian", "native_projected_hessian")):
                if (hessian_folder / source).exists():
                    attempt[stem] = _archive(hessian_folder / source, directory / f"{stem}_{attempt_index:02d}.gz")
            attempt["native_command"] = json.loads((hessian_folder / "command.json").read_text())
            frequencies = spectrum.frequencies_cm1
            negative = np.flatnonzero(frequencies < 0)
            retention = _coordination_check(current, original, coordination_indices)
            attempt.update(energy_eV=energy, native_hessian_wall_seconds=wall,
                frequencies_cm1=frequencies.tolist(), imaginary_count=int(len(negative)),
                external_modes_projected=spectrum.external_rank,
                hessian_antisymmetry_relative=spectrum.antisymmetry_relative,
                gradient_hessian_energy_difference_eV=abs(energy - energy_gradient), coordination=retention)
            if retention["checked"] and not retention["retained"]:
                attempt["status"] = "coordination_not_retained"
            elif force > force_threshold_ev_a:
                attempt["status"] = "forces_not_converged"
            elif spectrum.antisymmetry_relative > .05 or abs(energy - energy_gradient) > .01 or np.any(np.abs(frequencies) < 1e-5):
                attempt["status"] = "numerical_uncertainty"
            elif transition_state:
                in_window = (len(negative) == 1 and (ts_imaginary_window_cm1 is None or
                             ts_imaginary_window_cm1[0] <= frequencies[negative[0]] <= ts_imaginary_window_cm1[1]))
                attempt["status"] = ("verified_first_order_saddle" if in_window
                                     else "saddle_spectrum_rejected")
            elif len(negative):
                attempt["status"] = "numerical_uncertainty" if np.min(frequencies) >= -20 else "not_a_minimum"
            else:
                attempt["status"] = "verified_minimum"
            if attempt["status"] in ("verified_minimum", "verified_first_order_saddle"):
                records = thermochemistry_from_modes(current, energy, frequencies,
                    transition_state=transition_state, unpaired=unpaired, symmetry_number=symmetry_number,
                    temperatures=temperatures, solvent=solvent, solvation_state=solvation_state,
                    ts_imaginary_window_cm1=ts_imaginary_window_cm1,
                    source=str(spectrum_path),
                    source_sha256=attempt["spectrum_npz_sha256"])
                result.update(status=attempt["status"], accepted=True,
                              thermochemistry=[asdict(record) for record in records], accepted_attempt=attempt_index)
                write(directory / "accepted_geometry.xyz", current)
                result["accepted_geometry_sha256"] = hashlib.sha256((directory / "accepted_geometry.xyz").read_bytes()).hexdigest()
                break
            result["status"] = attempt["status"]
            atomic_json(directory / "result.json", result)
            if transition_state or attempt_index == max_repairs or attempt["status"] == "coordination_not_retained":
                break
            displaced = current.copy()
            if len(negative):
                mode = spectrum.cartesian_modes[negative[0]]
                displacement = 0.10 * (-1 if attempt_index % 2 else 1)
                displaced.positions += displacement * mode / np.max(np.linalg.norm(mode, axis=1))
                attempt["repair"] = {"mode_frequency_cm1": float(frequencies[negative[0]]), "max_displacement_A": displacement}
            else:
                attempt["repair"] = {"reason": "tighter native optimization without mode displacement"}
            repair_folder = scratch / f"repair_{attempt_index:02d}"
            attempt["phase"] = "repair_optimization"
            atomic_json(directory / "result.json", result)
            repair = optimize_geometry(displaced, repair_folder, charge=charge, unpaired=unpaired, threads=threads,
                executable=executable, max_cycles=200, convergence="extreme", solvent=solvent,
                solvation_state=solvation_state, timeout_seconds=_remaining(deadline_epoch, 900))
            attempt["repair"]["optimization"] = asdict(repair)
            if not repair.converged:
                result["status"] = "repair_not_converged"
                break
            current = read(repair_folder / "xtbopt.xyz")
            if not np.array_equal(current.numbers, original.numbers):
                raise ValueError("native repair changed atom identities or ordering")
    except (TimeoutError, subprocess.TimeoutExpired) as exc:
        result.update(status="deadline_exhausted", error=str(exc))
    except Exception as exc:
        result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    finally:
        if result["attempts"] and result["attempts"][-1]["status"] == "started":
            result["attempts"][-1].update(status=result["status"], error=result.get("error"))
        result.update(finished_utc=utc_now(), wall_seconds=time.monotonic() - started)
        atomic_json(directory / "result.json", result)
    return result
