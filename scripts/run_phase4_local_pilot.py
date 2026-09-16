"""A bounded, immutable three-job SCF recovery/basis diagnostic; no production labels.

The original 120 s quintet job was stopped one SCF energy tolerance from
convergence. This run keeps those tolerances and increases its time allowance.
Two subsequent same-geometry def2-SVP jobs test feasibility, not method accuracy.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import psutil
from ase.io import read
from ase.units import Bohr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.spin_provider import validate_spin_system

CATALYST = "Fe_bipyridine_pnnoh_iPr"
JOBS = (("sto3g_quintet_recovery", "sto-3g", 5, 240.0),
        ("def2svp_singlet", "def2-svp", 1, 450.0),
        ("def2svp_quintet", "def2-svp", 5, 450.0))
PROTOCOL = dict(method="pbe", scf_max_iterations=120, grid_radial=35,
                grid_spherical=110, scf_algorithm="soscf",
                soscf_start_convergence=1e-4, e_convergence=1e-8,
                d_convergence=1e-6, memory_mib=500, threads=2)


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def dump(path, value):
    target = Path(path)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(target)


def geometry_digest(symbols, positions):
    return hashlib.sha256(np.asarray(positions, dtype="<f8").tobytes() + " ".join(symbols).encode()).hexdigest()


def paths_for_run(output, scratch):
    output, scratch = Path(output).resolve(), Path(scratch).resolve()
    for path in (output, scratch):
        if path.exists():
            raise FileExistsError(f"Fresh path required; existing evidence preserved: {path}")
        if path == REPO / "data" / "phase3" or (REPO / "data" / "phase3") in path.parents:
            raise ValueError("Phase3 evidence is immutable")
    if output == scratch or output in scratch.parents or scratch in output.parents:
        raise ValueError("Output and scratch must be separate, non-nested directories")
    return output, scratch


def validate_response(response, request):
    """Reject geometry/protocol/state mismatches and non-finite worker results."""
    checks = {"symbols": request["symbols"], "charge": request["charge"],
              "multiplicities": request["multiplicities"], "method": request["method"],
              "basis": request["basis"], "coordinate_units": "bohr", "energy_units": "hartree",
              "gradient_units": "hartree/bohr", "scf_e_convergence": request["e_convergence"],
              "scf_d_convergence": request["d_convergence"]}
    checks.update({key: request[key] for key in ("grid_radial", "grid_spherical", "scf_algorithm",
                                               "scf_max_iterations", "soscf_start_convergence")})
    for key, value in checks.items():
        if response.get(key) != value:
            raise ValueError(f"Worker response differs from requested {key}")
    x = np.asarray(request["positions_bohr"], dtype=float)
    if response.get("geometry_sha256") != geometry_digest(request["symbols"], x):
        raise ValueError("Worker geometry hash mismatch")
    if not np.array_equal(np.asarray(response.get("positions_bohr")), x):
        raise ValueError("Worker coordinates mismatch")
    states = response.get("states", [])
    if len(states) != 1:
        raise ValueError("Exactly one worker state required")
    state = states[0]
    multiplicity = request["multiplicities"][0]
    gradient = np.asarray(state.get("gradient_hartree_bohr"), dtype=float)
    if (state.get("multiplicity") != multiplicity or state.get("scf_converged") is not True
            or state.get("gradient_kind") != "analytic"
            or state.get("reference") != ("rks" if multiplicity == 1 else "uks")):
        raise ValueError("Unconverged, wrong-reference, or wrong-state response")
    if gradient.shape != x.shape or not np.all(np.isfinite(gradient)):
        raise ValueError("Invalid analytic gradient")
    energy, s2 = float(state["energy_hartree"]), float(state["spin_squared"])
    if not np.isfinite(energy) or not np.isfinite(s2):
        raise ValueError("Nonfinite energy or spin")
    expected = (multiplicity * multiplicity - 1) / 4
    if state.get("expected_spin_squared") != expected:
        raise ValueError("Wrong expected spin")
    contamination = s2 - expected
    return dict(energy_hartree=energy, spin_squared=s2, expected_spin_squared=expected,
                spin_contamination=contamination, spin_quality_pass=abs(contamination) <= 0.1,
                max_gradient_hartree_bohr=float(np.linalg.norm(gradient, axis=1).max()),
                gradient_kind="analytic", scf_converged=True)


def archive_native(scratch, output):
    """Copy scientific text/arrays; retain/hash scratch caches without deleting them."""
    native = output / "native"
    native.mkdir(exist_ok=False)
    files = []
    for path in sorted(Path(scratch).rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(scratch)
        preserved = path.suffix.lower() in {".json", ".out", ".xyz", ".npy", ".npz"}
        item = dict(relative_path=rel.as_posix(), sha256=digest(path), bytes=path.stat().st_size,
                    scientific_copy=preserved, original_local_path=str(path))
        if preserved:
            destination = native / rel
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, destination)
            if digest(destination) != item["sha256"]:
                raise RuntimeError("Native evidence copy changed")
            item["published_relative_path"] = str(destination.relative_to(output)).replace("\\", "/")
        files.append(item)
    dump(output / "native_manifest.json", dict(files=files, scratch_files_deleted=False,
         note="Scientific text/arrays copied; PSIO binary intermediates retained locally with hashes."))
    return {"files": len(files), "scientific_files_copied": sum(f["scientific_copy"] for f in files),
            "manifest_sha256": digest(output / "native_manifest.json")}


def run_job(job, args, symbols, positions, output, scratch, deadline):
    name, basis, multiplicity, allowance = job
    folder, work = output / name, scratch / name
    folder.mkdir(exist_ok=False)
    work.mkdir(exist_ok=False)
    request = dict(**PROTOCOL, symbols=symbols, positions_bohr=positions.tolist(), charge=0,
                   multiplicities=[multiplicity], single_state=True, basis=basis,
                   workdir=str(work / "native"))
    validate_spin_system(symbols, positions, 0, [multiplicity], require_pair=False)
    dump(folder / "request.json", request)
    command = [str(args.python), "-m", "pincer_catmech.quantum.spin_provider", "--request",
               str(folder / "request.json"), "--result", str(folder / "response.json")]
    env = dict(os.environ)
    prefix = args.python.parent
    env.update(OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
               PYTHONUTF8="1", PYTHONIOENCODING="utf-8", PYTHONPATH=str(REPO / "src"))
    env["PATH"] = str(prefix / "Library" / "bin") + os.pathsep + str(prefix) + os.pathsep + env.get("PATH", "")
    started = time.monotonic()
    remaining = min(allowance, deadline - started)
    result = dict(name=name, basis=basis, method="pbe", multiplicity=multiplicity,
                  started_utc=utc(), status="not_started", requested_wall_seconds=allowance,
                  enforced_wall_seconds=max(0, remaining), configured_memory_mib=500,
                  rss_limit_mib=args.rss_limit_mib, threads=2, peak_rss_mib=0.0,
                  accepted_chemical_label=False, geometry_sha256=geometry_digest(symbols, positions))
    process = None
    try:
        if remaining < 30:
            result.update(status="not_started_budget", reason="Under 30 seconds remain")
            return result
        if psutil.virtual_memory().available < (args.rss_limit_mib + 512) * 1024**2:
            result.update(status="not_started_memory", reason="Available memory below RSS limit plus 512 MiB reserve")
            return result
        with (folder / "worker.out").open("w", encoding="utf-8") as log:
            process = subprocess.Popen(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
            measured = psutil.Process(process.pid)
            result["status"] = "running"
            while process.poll() is None:
                try:
                    rss = measured.memory_info().rss + sum(c.memory_info().rss for c in measured.children(recursive=True))
                    result["peak_rss_mib"] = max(result["peak_rss_mib"], rss / 1024**2)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    rss = 0
                if time.monotonic() - started >= remaining:
                    result.update(status="timeout", reason="Enforced per-job or total wall-time limit")
                    break
                if rss > args.rss_limit_mib * 1024**2 or psutil.virtual_memory().available < 256 * 1024**2:
                    result.update(status="memory_limit", reason="RSS limit or host free-memory reserve reached")
                    break
                time.sleep(0.1)
            if process.poll() is None:
                for child in measured.children(recursive=True):
                    child.kill()
                process.kill()
            result["returncode"] = process.wait()
        if result["status"] == "running":
            if result["returncode"] != 0 or not (folder / "response.json").is_file():
                result.update(status="worker_failed", reason="Nonzero exit or missing response")
            else:
                response = json.loads((folder / "response.json").read_text(encoding="utf-8"))
                result.update(validate_response(response, request), status="converged")
    except Exception as error:
        result.update(status="driver_failed", reason=f"{type(error).__name__}: {error}")
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
        result.update(finished_utc=utc(), worker_wall_seconds=time.monotonic() - started)
        result["native_evidence"] = archive_native(work, folder)
        result["command"] = command
        result["request_sha256"] = digest(folder / "request.json")
        if (folder / "worker.out").is_file():
            result["worker_stdout_sha256"] = digest(folder / "worker.out")
        if (folder / "response.json").is_file():
            result["response_sha256"] = digest(folder / "response.json")
        dump(folder / "result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True, type=Path, help="Existing Psi4 environment Python")
    parser.add_argument("--output", type=Path, default=REPO / "data/phase4/local_pilot_001")
    parser.add_argument("--scratch", type=Path, default=REPO.parent / "phase4/local_pilot_001")
    parser.add_argument("--wall-seconds", type=float, default=1200)
    parser.add_argument("--rss-limit-mib", type=int, default=1536)
    args = parser.parse_args()
    if not np.isfinite(args.wall_seconds) or not 60 <= args.wall_seconds <= 1200:
        parser.error("Total wall budget must be finite, 60..1200 seconds")
    if not 600 <= args.rss_limit_mib <= 1536:
        parser.error("RSS limit must be 600..1536 MiB")
    args.python = args.python.resolve()
    if not args.python.is_file():
        raise FileNotFoundError(args.python)
    output, scratch = paths_for_run(args.output, args.scratch)
    source = REPO / "data/structures" / CATALYST / "active/best_found.xyz"
    metadata_path = source.with_name("metadata.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    atoms = read(source)
    identity = catalyst_identity_screen(atoms, metadata)
    if metadata.get("charge") != 0 or not identity.get("retained"):
        raise ValueError("Input charge or mapped chemical identity failed")
    symbols, positions = atoms.get_chemical_symbols(), atoms.positions / Bohr
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    sources = output / "sources"
    sources.mkdir()
    files = [Path(__file__), REPO / "src/pincer_catmech/quantum/spin_provider.py", source, metadata_path]
    source_records = []
    for path in files:
        target = sources / path.name
        shutil.copyfile(path, target)
        source_records.append(dict(original=str(path), copy=str(target.relative_to(output)), sha256=digest(target)))
    engine_files = [args.python] + sorted((args.python.parent / "conda-meta").glob("psi4-*.json"))
    engine_files += sorted((args.python.parent / "Lib/site-packages/psi4").glob("core*.pyd"))
    summary = dict(schema="phase4_local_spin_pilot_v1", started_utc=utc(), status="running", catalyst_id=CATALYST,
                   source_files=source_records, engine_fingerprint={str(p): digest(p) for p in engine_files},
                   identity=identity, jobs=[], total_wall_budget_seconds=args.wall_seconds,
                   physical_workers=1, threads=2, baseline_temperature_K=383.15,
                   electronic_model_temperature_K=None, solvent=None, geometry_optimized_by_this_run=False,
                   scope="Gas-phase vertical SCF and basis feasibility diagnostic at a prior xTB geometry.",
                   production_spin_gap_validated=False, accepted_MECP_count=0)
    started = time.monotonic()
    deadline = started + args.wall_seconds
    dump(output / "summary.json", summary)
    try:
        for job in JOBS:
            print(f"Starting {job[0]} with up to {job[3]:.0f} seconds", flush=True)
            result = run_job(job, args, symbols, positions, output, scratch, deadline)
            summary["jobs"].append(result)
            dump(output / "summary.json", summary)
            print(json.dumps({key: result.get(key) for key in ("name", "status", "energy_hartree", "spin_squared", "worker_wall_seconds")}), flush=True)
        summary["status"] = "completed_bounded_diagnostic"
    finally:
        summary.update(finished_utc=utc(), elapsed_seconds=time.monotonic() - started,
                       scf_converged_count=sum(j["status"] == "converged" for j in summary["jobs"]),
                       spin_quality_pass_count=sum(j.get("spin_quality_pass", False) for j in summary["jobs"]))
        dump(output / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
