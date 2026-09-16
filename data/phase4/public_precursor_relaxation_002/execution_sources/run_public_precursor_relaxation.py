"""Bounded native xTB preprocessing of two complete public Mn precursor crystals.

Two optimizations run first; complete Hessians follow only within remaining
budgets. No retry, ground-spin claim, activation chemistry, or DFT validation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import psutil
from ase.data import covalent_radii
from ase.io import read, write

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from pincer_catmech.quantum.xtb_backend import (
    _execute, atomic_json, executable_path, optimize_geometry, parse_energy,
    parse_gradient, utc_now, validate_electrons,
)
from pincer_catmech.quantum.thermochemistry import native_hessian

SOURCE = REPO / "data/phase4/public_structure/import_v002"
PER_STRUCTURE_SECONDS = 420.0
OPT_SECONDS = 210.0
FORCE_THRESHOLD = .03


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fresh_paths(output, scratch):
    output, scratch = Path(output).resolve(), Path(scratch).resolve()
    if output.exists() or scratch.exists():
        raise FileExistsError("Fresh output and scratch required; prior evidence remains unchanged")
    if output == scratch or output in scratch.parents or scratch in output.parents:
        raise ValueError("Output and scratch must be separate and non-nested")
    for base in (REPO / "data/phase3", SOURCE):
        if any(path == base or base in path.parents for path in (output, scratch)):
            raise ValueError("Historical evidence and source crystal paths are immutable")
    return output, scratch


def graph_identity(atoms, source, metadata):
    """CIF-derived atom mapping and all bonds, not generated-catalog chemistry."""
    if (len(atoms) != 62 or not np.array_equal(atoms.numbers, source.numbers)
            or not np.all(np.isfinite(atoms.positions))):
        return dict(retained=False, reason="Atom count, ordered elements or finite geometry mismatch")
    expected = {tuple(sorted(pair)) for pair in metadata["all_bond_indices"]}
    if any(a == b or a < 0 or b >= len(atoms) for a, b in expected):
        raise ValueError("Invalid mapped source bond")
    distances = atoms.get_all_distances(mic=False)
    actual = {(a, b) for a in range(len(atoms)) for b in range(a + 1, len(atoms))
              if distances[a, b] < 1.20 * (covalent_radii[atoms.numbers[a]] + covalent_radii[atoms.numbers[b]])}
    m, n, h = metadata["metal_index"], metadata["proton_site_index"], metadata["ligand_proton_index"]
    essential = [tuple(sorted((m, i))) for i in metadata["donor_indices"] + metadata["halide_indices"]]
    essential += [tuple(sorted((m, c))) for c, o in metadata["carbonyl_indices"]]
    essential += [tuple(sorted(pair)) for pair in metadata["carbonyl_indices"]] + [tuple(sorted((n, h)))]
    if (atoms[m].symbol != "Mn" or atoms[n].symbol != "N" or atoms[h].symbol != "H"
            or len(metadata["carbonyl_indices"]) != 2 or len(metadata["halide_indices"]) != 1
            or sorted(atoms[i].symbol for i in metadata["donor_indices"]) != ["N", "P", "P"]):
        raise ValueError("Source is not the expected NH/P2N/Br/(CO)2 precursor")
    return dict(retained=actual == expected and all(pair in actual for pair in essential),
                lost_bonds=[list(x) for x in sorted(expected - actual)], new_bonds=[list(x) for x in sorted(actual - expected)],
                essential_bonds_retained=all(pair in actual for pair in essential),
                NH_distance_angstrom=float(distances[n, h]),
                Mn_Br_distance_angstrom=float(distances[m, metadata["halide_indices"][0]]),
                Mn_donor_distances_angstrom=[float(distances[m, i]) for i in metadata["donor_indices"]],
                carbonyl_distances_angstrom=[dict(Mn_C=float(distances[m, c]), C_O=float(distances[c, o])) for c, o in metadata["carbonyl_indices"]],
                criterion="Exact mapped source bond graph under 1.20 times ASE covalent-radii sum; geometric screen, not bond-order or solution-identity proof")


def remaining_budget(record, deadline):
    return max(0., min(PER_STRUCTURE_SECONDS - record["quantum_wall_seconds"], deadline - time.time()))


def validate_spectrum(spectrum, atom_count):
    frequencies = np.asarray(spectrum.frequencies_cm1, dtype=float)
    if frequencies.shape != (3 * atom_count - 6,) or not np.all(np.isfinite(frequencies)):
        raise ValueError("Full nonlinear-molecule internal mode cardinality required")
    return dict(full_internal_mode_count=len(frequencies), imaginary_count=int(np.count_nonzero(frequencies < 0)),
                zero_mode_count=int(np.count_nonzero(frequencies == 0)), lowest_frequency_cm1=float(frequencies.min()),
                hessian_antisymmetry_relative=float(spectrum.antisymmetry_relative),
                hessian_minimum_pass=bool(np.all(frequencies > 0) and np.isfinite(spectrum.antisymmetry_relative)
                                         and spectrum.antisymmetry_relative <= .05))


def preserve_native(work, folder):
    target = folder / "native"
    shutil.copytree(work, target)
    rows = []
    for path in sorted(work.rglob("*")):
        if path.is_file():
            relative = path.relative_to(work)
            checksum = digest(path)
            if digest(target / relative) != checksum:
                raise RuntimeError("Native evidence copy mismatch")
            rows.append(dict(path=(Path("native") / relative).as_posix(), sha256=checksum, bytes=path.stat().st_size))
    atomic_json(folder / "native_manifest.json", dict(files=rows, all_native_files_copied=True, scratch_deleted=False))
    return dict(file_count=len(rows), manifest_sha256=digest(folder / "native_manifest.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb", default=os.environ.get("PINCER_XTB"))
    parser.add_argument("--output", type=Path, default=REPO / "data/phase4/public_precursor_relaxation_001")
    parser.add_argument("--scratch", type=Path, default=REPO.parent / "phase4/public_precursor_relaxation_001")
    parser.add_argument("--wall-seconds", type=float, default=900)
    args = parser.parse_args()
    if not np.isfinite(args.wall_seconds) or not 60 <= args.wall_seconds <= 900:
        parser.error("Use finite 60..900 second total quantum budget")
    executable = executable_path(args.xtb)
    output, scratch = fresh_paths(args.output, args.scratch)
    manifest = json.loads((SOURCE / "file_manifest.json").read_text(encoding="utf-8"))
    seeds = []
    for label in ("Mn1", "Mn2"):
        xyz, meta_path = SOURCE / label / "crystal_candidate.xyz", SOURCE / label / "metadata.json"
        for path in (xyz, meta_path):
            if digest(path) != manifest[path.relative_to(SOURCE).as_posix()]:
                raise ValueError("Imported source geometry/metadata no longer match frozen manifest")
        meta, atoms = json.loads(meta_path.read_text(encoding="utf-8")), read(xyz)
        if meta["charge"] != 0 or meta["unpaired_electrons"] != 0:
            raise ValueError("Explicit q0/uhf0 candidate source required")
        validate_electrons(atoms, 0, 0)
        check = graph_identity(atoms, atoms, meta)
        if not check["retained"]:
            raise ValueError("Original crystal candidate fails its own full mapped graph screen")
        seeds.append((label, atoms, meta, xyz, meta_path, check))
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    sources = output / "execution_sources"
    sources.mkdir()
    source_paths = [Path(__file__), REPO / "src/pincer_catmech/quantum/xtb_backend.py",
                    REPO / "src/pincer_catmech/quantum/thermochemistry.py",
                    REPO / "src/pincer_catmech/kinetics/neb_ts_search.py"]
    source_hashes = {}
    for path in source_paths:
        shutil.copyfile(path, sources / path.name)
        source_hashes[path.relative_to(REPO).as_posix()] = digest(path)
    summary = dict(schema="public_precursor_relaxation_v1", status="running", started_utc=utc_now(),
                   method="GFN2-xTB", solvent="toluene", solvation_state="gsolv", threads=2, physical_workers=1,
                   charge=0, unpaired_electrons=0, electronic_smearing_temperature_K=300.,
                   reaction_temperature_context_K=383.15, thermal_free_energy_computed=False,
                   configured_memory_mib=None, memory_limit_enforced=False,
                   memory_policy="Require at least 1 GiB available host memory before each native calculation; native xTB has no imposed allocation here",
                   total_quantum_wall_limit_seconds=args.wall_seconds, per_structure_quantum_limit_seconds=420,
                   optimization_limit_seconds=210, force_threshold_eV_A=FORCE_THRESHOLD,
                   hessian_minimum_remaining_seconds=120, retries=0, structures=[], execution_source_sha256=source_hashes,
                   xtb_executable=str(executable), xtb_executable_sha256=digest(executable),
                   source_cif_sha256=digest(SOURCE / "original.cif"),
                   ground_spin_determined=False, validated_solution_mechanism=False,
                   new_transition_state_count=0, new_MECP_count=0,
                   scope="xTB preparation/identity diagnostic on two independent crystal conformers of the same precursor; no ground-spin, active-species, DFT or experimental validation")
    start = time.monotonic()
    deadline = time.time() + args.wall_seconds
    atomic_json(output / "summary.json", summary)
    optimized = {}
    try:
        # Both optimizations precede Hessians, so one Hessian cannot starve the second seed.
        for label, atoms, metadata, xyz, meta_path, check in seeds:
            folder, work = output / label, scratch / label
            folder.mkdir()
            work.mkdir()
            shutil.copyfile(xyz, folder / "crystal_input.xyz")
            shutil.copyfile(meta_path, folder / "source_metadata.json")
            record = dict(id=label, status="not_started", source_xyz_sha256=digest(xyz), source_metadata_sha256=digest(meta_path),
                          input_identity=check, optimization_converged=False, minimum_certified=False,
                          accepted_chemical_label=False, quantum_wall_seconds=0., native_calculations_started=0,
                          failure_reasons=[], hessian_status="not_started")
            summary["structures"].append(record)
            try:
                allowed = min(OPT_SECONDS, remaining_budget(record, deadline))
                if allowed < 30 or psutil.virtual_memory().available < 1024**3:
                    record.update(status="not_started_resource", reason="Insufficient remaining time or 1 GiB memory reserve")
                    continue
                print(f"Optimizing {label}, at most {allowed:.1f} s", flush=True)
                quantum_start = time.monotonic()
                record["native_calculations_started"] += 1
                result = optimize_geometry(atoms, work / "optimization", charge=0, unpaired=0, threads=2,
                    executable=executable, max_cycles=500, timeout_seconds=allowed, convergence="extreme",
                    solvent="toluene", solvation_state="gsolv")
                record["quantum_wall_seconds"] += time.monotonic() - quantum_start
                record["optimization"] = asdict(result)
                record["optimization_converged"] = result.converged
                if not result.converged:
                    record["status"] = "optimization_failed"
                    continue
                candidate = read(work / "optimization/xtbopt.xyz")
                identity = graph_identity(candidate, atoms, metadata)
                record["relaxed_identity"] = identity
                write(folder / "relaxed_candidate.xyz", candidate, format="xyz")
                record["relaxed_candidate_sha256"] = digest(folder / "relaxed_candidate.xyz")
                record["status"] = "optimized_unverified_minimum" if identity["retained"] else "identity_changed"
                if identity["retained"]:
                    optimized[label] = candidate
            except Exception as error:
                record["status"] = "optimization_driver_failed"
                record["failure_reasons"].append(f"{type(error).__name__}: {error}")
            finally:
                atomic_json(folder / "result.json", record)
                atomic_json(output / "summary.json", summary)
        for record in summary["structures"]:
            label = record["id"]
            if label not in optimized:
                continue
            atoms = optimized[label]
            work, folder = scratch / label, output / label
            try:
                allowed = min(30., remaining_budget(record, deadline))
                if allowed < 10 or psutil.virtual_memory().available < 1024**3:
                    record["hessian_status"] = "not_started_resource"
                    continue
                quantum_start = time.monotonic()
                record["native_calculations_started"] += 1
                try:
                    code, text, wall, argv = _execute(atoms, work / "fresh_gradient", executable, 0, 0, 2,
                        ["--grad"], allowed, electronic_temperature=300., solvent="toluene", solvation_state="gsolv")
                    if code != 0 or "normal termination of xtb" not in text.lower():
                        raise ValueError("Fresh gradient calculation did not terminate normally")
                    forces = parse_gradient(work / "fresh_gradient/gradient", len(atoms))
                    record["fresh_max_force_eV_A"] = float(np.linalg.norm(forces, axis=1).max())
                    record["electronic_energy_eV"] = parse_energy(text)
                finally:
                    record["quantum_wall_seconds"] += time.monotonic() - quantum_start
                if record["fresh_max_force_eV_A"] > FORCE_THRESHOLD:
                    record["hessian_status"] = "not_started_force_failed"
                    continue
                allowed = remaining_budget(record, deadline)
                if allowed < 120 or psutil.virtual_memory().available < 1024**3:
                    record["hessian_status"] = "not_started_resource"
                    continue
                print(f"Full Hessian {label}, at most {allowed:.1f} s", flush=True)
                quantum_start = time.monotonic()
                record["native_calculations_started"] += 1
                try:
                    spectrum, energy, wall = native_hessian(atoms, work / "hessian", charge=0, unpaired=0, threads=2,
                        executable=executable, deadline_epoch=deadline, timeout_seconds=allowed,
                        solvent="toluene", solvation_state="gsolv")
                    spectrum.save(folder / "full_hessian.npz", atoms)
                    record.update(validate_spectrum(spectrum, len(atoms)), hessian_status="completed",
                                  hessian_energy_eV=energy, frequencies_cm1=spectrum.frequencies_cm1.tolist(),
                                  full_hessian_sha256=digest(folder / "full_hessian.npz"))
                    record["minimum_certified"] = record["hessian_minimum_pass"]
                    if record["minimum_certified"]:
                        record["status"] = "verified_model_minimum"
                finally:
                    record["quantum_wall_seconds"] += time.monotonic() - quantum_start
            except Exception as error:
                record["hessian_status"] = "failed"
                record["failure_reasons"].append(f"{type(error).__name__}: {error}")
            finally:
                atomic_json(folder / "result.json", record)
                atomic_json(output / "summary.json", summary)
        summary["status"] = "completed_bounded_preprocessing"
    finally:
        for record in summary["structures"]:
            label = record["id"]
            record["native_evidence"] = preserve_native(scratch / label, output / label)
            atomic_json(output / label / "result.json", record)
        summary.update(finished_utc=utc_now(), total_elapsed_seconds=time.monotonic() - start,
                       optimizations_converged=sum(r["optimization_converged"] for r in summary["structures"]),
                       mapped_structures_retained=sum(r.get("relaxed_identity", {}).get("retained", False) for r in summary["structures"]),
                       verified_model_minima=sum(r["minimum_certified"] for r in summary["structures"]),
                       native_calculations_started=sum(r["native_calculations_started"] for r in summary["structures"]))
        atomic_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
