"""Real, bounded condensation CI-NEB with mechanism-specific certification.

This module reuses numerical Hessian/dimer tools without imposing the alcohol
dehydrogenation hydride criterion. A converged band is not a certified saddle.
"""
from dataclasses import asdict
from pathlib import Path
import hashlib
import time

import numpy as np
from ase.data import covalent_radii
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE

from pincer_catmech.kinetics.neb_ts_search import refine_dimer
from pincer_catmech.quantum.thermochemistry import native_hessian, thermochemistry_from_modes
from pincer_catmech.quantum.xtb_backend import XTBCalculator, atomic_json, optimize_geometry, utc_now


def observed_bonds(atoms):
    return {tuple((j, i)) for i in range(len(atoms)) for j in range(i)
            if atoms.get_distance(i, j) < 1.28*(covalent_radii[atoms.numbers[i]]+covalent_radii[atoms.numbers[j]])}


def identity(atoms, expected_bonds, *, cn, cn_interval):
    expected = {tuple(sorted(pair)) for pair in expected_bonds}
    actual = observed_bonds(atoms)
    distance = float(atoms.get_distance(*cn))
    return {"retained": expected == actual and cn_interval[0] <= distance <= cn_interval[1],
            "lost_bonds": sorted(expected-actual), "new_bonds": sorted(actual-expected),
            "C_N_distance_angstrom": distance, "C_N_allowed_interval": list(cn_interval),
            "criterion": "1.28 times ASE covalent radius sum plus specified C-N distance interval"}


def coordinate_gradient(atoms, terms):
    gradient = np.zeros_like(atoms.positions)
    for i, j, coefficient in terms:
        if not (0 <= i < len(atoms) and 0 <= j < len(atoms)) or i == j or not np.isfinite(coefficient):
            raise ValueError("Invalid mapped reaction-coordinate pair or coefficient")
        difference = atoms.positions[i]-atoms.positions[j]
        norm = np.linalg.norm(difference)
        if norm < 1e-8:
            raise ValueError("Coincident atoms in reaction coordinate")
        gradient[i] += coefficient*difference/norm
        gradient[j] -= coefficient*difference/norm
    return gradient


def condensation_mode(atoms, mode, coordinates):
    mode = np.asarray(mode, dtype=float)
    if mode.shape != atoms.positions.shape or not np.all(np.isfinite(mode)) or np.linalg.norm(mode) < 1e-12 or not coordinates:
        raise ValueError("A finite nonzero full Cartesian mode and reaction coordinates are required")
    unit = mode/np.linalg.norm(mode)
    overlaps = {}
    total = np.zeros_like(mode)
    for name, spec in coordinates.items():
        gradient = coordinate_gradient(atoms, spec["terms"])
        if np.linalg.norm(gradient) < 1e-12:
            raise ValueError("Reaction coordinate has a vanishing Cartesian gradient")
        gradient /= np.linalg.norm(gradient)
        overlaps[name] = float(np.sum(unit*gradient))
        total += gradient
    if np.linalg.norm(total) < 1e-12:
        raise ValueError("Reaction coordinate directions cancel")
    sign = 1 if np.sum(unit*total) >= 0 else -1
    aligned = {name: sign*value for name, value in overlaps.items()}
    combined = float(abs(np.sum(unit*total))/np.linalg.norm(total))
    accepted = all(aligned[name] >= spec["minimum_overlap"] for name, spec in coordinates.items()) and combined >= 0.25
    return {"accepted": accepted, "signed_overlaps": overlaps, "oriented_overlaps": aligned,
            "combined_overlap": combined, "criterion": "same direction for all explicitly declared proton-transfer/bond-change coordinates; no hydride requirement"}


def run_condensation_path(reactant, product, *, output_dir, scratch_dir, executable,
                          charge, endpoint_bonds, reaction_coordinates, cn,
                          deadline_epoch, neb_steps=180, endpoint_cycles=400,
                          neb_fmax=0.07, ts_fmax=0.03, model_label,
                          restart_band=None, neb_maxstep=0.10, certify_endpoints=False):
    """Nine total images, actual ALPB energies, dimer, full Hessian, two descents.

    Outputs have accepted=False unless all chemical/numerical gates pass.
    Electronic differences and modeled qRRHO quantities remain distinguished.
    """
    if not np.isfinite(neb_maxstep) or neb_maxstep <= 0 or neb_steps < 1:
        raise ValueError("Positive finite NEB maximum step and positive step budget required")
    directory, scratch = Path(output_dir), Path(scratch_dir)
    directory.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    if not np.array_equal(reactant.numbers, product.numbers):
        raise ValueError("Both endpoints must retain exactly the same mapped atom order")
    start = time.time()
    result = {"status": "running", "accepted": False, "started_utc": utc_now(),
              "model": model_label, "charge": charge, "unpaired_electrons": 0,
              "method": "GFN2-xTB", "solvent": "toluene", "solvation_state": "gsolv",
              "threads": 1, "n_total_images": 9, "n_internal_images": 7,
              "atom_mapping": list(range(len(reactant))), "endpoints": [], "stages": [],
              "reaction_coordinates": reaction_coordinates, "failure_reasons": [],
              "free_energy_barrier_computed": False}
    restart = None
    if restart_band is not None:
        source = Path(restart_band).resolve()
        restart = read(source, index=":")
        if len(restart) != 9 or any(not np.array_equal(frame.numbers, reactant.numbers) for frame in restart):
            raise ValueError("Restart requires nine complete frames in identical atom order")
        if not np.allclose(restart[0].positions, reactant.positions, rtol=0, atol=1e-8) or not np.allclose(restart[-1].positions, product.positions, rtol=0, atol=1e-8):
            raise ValueError("Restart endpoint coordinates must match the supplied endpoints")
        result["restart"] = {"source": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                             "all_forces_recomputed": True, "new_interpolation": False}
    atomic_json(directory/"protocol.json", {**result, "deadline_epoch": deadline_epoch,
        "neb_fmax_eV_A": neb_fmax, "ts_fmax_eV_A": ts_fmax, "neb_steps": neb_steps,
        "neb_maxstep_angstrom": neb_maxstep,
        "full_endpoint_hessians_required": certify_endpoints,
        "endpoint_cycles": endpoint_cycles, "imaginary_frequency_count_required": 1,
        "imaginary_frequency_window": None, "hydride_mode_required": False,
        "endpoint_bonds": endpoint_bonds, "C_N_indices": cn})

    def check():
        if time.time() >= deadline_epoch:
            raise TimeoutError("Condensation search wall-time allowance exhausted")

    def emit(stage, **values):
        entry = {"stage": stage, "utc": utc_now(), "elapsed_seconds": time.time()-start, **values}
        result["stages"].append(entry)
        atomic_json(directory/"progress.json", entry)
        atomic_json(directory/"result.json", result)
        print(entry, flush=True)

    def calculator(label):
        return XTBCalculator(scratch/label, executable=executable, charge=charge,
            unpaired=0, threads=1, solvent="toluene", solvation_state="gsolv",
            deadline=deadline_epoch, timeout_seconds=60)

    def optimize(atoms, label, cycles):
        check()
        outcome = optimize_geometry(atoms, scratch/label, charge=charge, unpaired=0,
            threads=1, executable=executable, solvent="toluene", solvation_state="gsolv",
            max_cycles=cycles, timeout_seconds=min(300, deadline_epoch-time.time()),
            convergence="extreme" if certify_endpoints else "tight")
        atomic_json(directory/(label+"_optimization.json"), asdict(outcome))
        if not outcome.converged:
            raise ValueError(label+" native geometry optimization did not converge")
        optimized = read(scratch/label/"xtbopt.xyz")
        optimized.calc = calculator(label+"_verify")
        force = float(np.linalg.norm(optimized.get_forces(), axis=1).max())
        if force > ts_fmax:
            raise ValueError(f"{label} has excessive unmodified force: {force}")
        write(directory/(label+".xyz"), optimized)
        return optimized, outcome

    try:
        endpoints = []
        for i, atoms in enumerate((reactant, product)):
            write(directory/f"endpoint_{i}_input.xyz", atoms)
            if restart is None:
                optimized, native = optimize(atoms, f"endpoint_{i}", endpoint_cycles)
                native_record = asdict(native)
            else:
                optimized = atoms.copy()
                optimized.calc = calculator(f"endpoint_{i}_restart_verify")
                force = float(np.linalg.norm(optimized.get_forces(), axis=1).max())
                if force > ts_fmax:
                    raise ValueError("Restart endpoint fails fresh unmodified force validation")
                native_record = {"status": "previous_optimized_coordinates_revalidated", "max_force_eV_A": force}
                write(directory/f"endpoint_{i}.xyz", optimized)
            basin = identity(optimized, endpoint_bonds[i], cn=cn, cn_interval=(1.35, 1.70) if i == 0 else (1.15, 1.42))
            result["endpoints"].append({"index": i, "energy_eV": float(optimized.get_potential_energy()),
                                       "identity": basin, "native": native_record})
            emit("endpoint", index=i, identity=basin)
            if not basin["retained"]:
                raise ValueError(f"Endpoint {i} optimized into another chemical basin")
            if certify_endpoints:
                spectrum, endpoint_energy, wall = native_hessian(optimized, scratch/f"endpoint_{i}_hessian",
                    charge=charge, executable=executable, threads=1, deadline_epoch=deadline_epoch,
                    solvent="toluene", solvation_state="gsolv")
                spectrum_path = directory/f"endpoint_{i}_hessian.npz"
                spectrum.save(spectrum_path, optimized)
                result["endpoints"][-1].update(frequencies_cm1=spectrum.frequencies_cm1.tolist(),
                    hessian_antisymmetry_relative=spectrum.antisymmetry_relative,
                    native_hessian_wall_seconds=wall, hessian_energy_eV=endpoint_energy,
                    hessian_sha256=hashlib.sha256(spectrum_path.read_bytes()).hexdigest())
                if np.any(spectrum.frequencies_cm1 <= 0) or spectrum.antisymmetry_relative > 0.05:
                    raise ValueError(f"Endpoint {i} full Hessian did not certify a minimum")
                result["endpoints"][-1]["stationarity_status"] = "verified_minimum"
                result["endpoints"][-1]["thermochemistry"] = [asdict(row) for row in thermochemistry_from_modes(
                    optimized, endpoint_energy, spectrum.frequencies_cm1, source=str(spectrum_path),
                    source_sha256=hashlib.sha256(spectrum_path.read_bytes()).hexdigest(),
                    solvent="toluene", solvation_state="gsolv")]
                emit("endpoint_hessian", index=i, minimum_frequency_cm1=float(spectrum.frequencies_cm1.min()))
            endpoints.append(optimized)
        images = [endpoints[0]]+([endpoints[0].copy() for _ in range(7)] if restart is None else [frame.copy() for frame in restart[1:-1]])+[endpoints[1]]
        for i, atoms in enumerate(images[1:-1], 1):
            atoms.calc = calculator(f"neb_{i:02d}")
        neb = NEB(images, climb=restart is not None, parallel=False, method="improvedtangent",
                  remove_rotation_and_translation=True, allow_shared_calculator=False)
        if restart is None:
            neb.interpolate(method="idpp", apply_constraint=False)
        write(directory/"initial_band.xyz", images)
        with FIRE(neb, logfile=str(directory/"neb.log"), trajectory=str(directory/"neb.traj"), maxstep=neb_maxstep) as optimizer:
            def progress():
                check()
                energies = [float(image.get_potential_energy()) for image in images]
                emit("neb", step=optimizer.nsteps, image_energies_eV=energies)
                write(directory/"current_band.xyz", images)
            optimizer.attach(progress, interval=5)
            if restart is None:
                optimizer.run(fmax=neb_fmax, steps=min(30, neb_steps//4))
            neb.climb = True
            converged = bool(optimizer.run(fmax=neb_fmax, steps=neb_steps-optimizer.nsteps))
        energies = [float(image.get_potential_energy()) for image in images]
        result.update(neb_converged=converged, image_energies_eV=energies)
        write(directory/"final_band.xyz", images)
        if not converged:
            raise ValueError("CI-NEB did not converge within the bounded step allowance")
        peak = int(np.argmax(energies[1:-1]))+1
        candidate = images[peak].copy()
        tangent = images[peak+1].positions-images[peak-1].positions
        candidate.calc = calculator("dimer")
        refined = refine_dimer(candidate, tangent, directory/"dimer", fmax=ts_fmax, steps=160)
        write(directory/"saddle_candidate.xyz", candidate)
        emit("dimer", converged=refined, true_fmax_eV_A=float(np.linalg.norm(candidate.get_forces(), axis=1).max()))
        if not refined:
            raise ValueError("Dimer refinement did not meet the true force criterion")
        spectrum, energy, wall = native_hessian(candidate, scratch/"native_hessian", charge=charge,
            executable=executable, threads=1, deadline_epoch=deadline_epoch, solvent="toluene", solvation_state="gsolv")
        spectrum.save(directory/"full_hessian.npz", candidate)
        if energy <= max(energies[0], energies[-1]):
            raise ValueError("Candidate saddle is not above both optimized endpoint clusters")
        imaginary = np.flatnonzero(spectrum.frequencies_cm1 < 0)
        result.update(frequencies_cm1=spectrum.frequencies_cm1.tolist(), imaginary_count=int(len(imaginary)),
            hessian_antisymmetry_relative=spectrum.antisymmetry_relative, hessian_wall_seconds=wall)
        emit("hessian", imaginary_count=int(len(imaginary)), frequencies_cm1=spectrum.frequencies_cm1.tolist())
        if len(imaginary) != 1 or spectrum.antisymmetry_relative > 0.05:
            raise ValueError("Full internal Hessian is not a numerically consistent first-order saddle")
        mode = spectrum.cartesian_modes[imaginary[0]]
        diagnostic = condensation_mode(candidate, mode, reaction_coordinates)
        result["reaction_mode"] = diagnostic
        if not diagnostic["accepted"]:
            raise ValueError("Imaginary mode does not describe the specified condensation/proton transfers")
        branches = []
        for sign, label in ((-1, "minus"), (1, "plus")):
            branch = candidate.copy()
            branch.positions += sign*0.20*mode/np.max(np.linalg.norm(mode, axis=1))
            relaxed, _ = optimize(branch, "descent_"+label, 500)
            checks = [identity(relaxed, bonds, cn=cn, cn_interval=(1.35, 1.70) if i == 0 else (1.15, 1.42))
                      for i, bonds in enumerate(endpoint_bonds)]
            matches = [i for i, basin in enumerate(checks) if basin["retained"]]
            branches.append({"sign": sign, "energy_eV": float(relaxed.get_potential_energy()),
                             "endpoint_matches": matches, "identity_checks": checks})
        result["descent_branches"] = branches
        if sorted(branch["endpoint_matches"] for branch in branches) != [[0], [1]] or any(branch["energy_eV"] >= energy-1e-5 for branch in branches):
            raise ValueError("Downhill displaced-mode relaxations did not connect opposite correct basins")
        source = directory/"full_hessian.npz"
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        result["thermochemistry"] = [asdict(row) for row in thermochemistry_from_modes(candidate,
            energy, spectrum.frequencies_cm1, transition_state=True, source=str(source),
            source_sha256=source_hash, solvent="toluene", solvation_state="gsolv",
            ts_imaginary_window_cm1=None)]
        result.update(status="verified_first_order_saddle", accepted=True, electronic_energy_eV=energy,
            forward_cluster_electronic_barrier_eV=energy-energies[0],
            reverse_cluster_electronic_barrier_eV=energy-energies[-1],
            source_sha256=source_hash, source_path=str(source), connectivity_verified=True)
        write(directory/"accepted_ts.xyz", candidate)
    except TimeoutError as error:
        result.update(status="deadline_exhausted")
        result["failure_reasons"].append(f"{type(error).__name__}: {error}")
    except Exception as error:
        result.update(status="failed")
        result["failure_reasons"].append(f"{type(error).__name__}: {error}")
    finally:
        result.update(finished_utc=utc_now(), elapsed_seconds=time.time()-start)
        atomic_json(directory/"result.json", result)
    return result
