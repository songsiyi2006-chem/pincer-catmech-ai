"""Bounded, prospective pincer MECP search from an eligible completed spin matrix.

No vertical gap or its sign certifies a crossing. With no eligible matched pair,
write an explicit readiness receipt and make zero electronic-structure calls.
The default search and optional curvature check share eight pair evaluations.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path, PureWindowsPath
import subprocess
import sys
import time
import uuid

import numpy as np
import psutil
from ase import Atoms
from ase.io import read
from ase.units import Bohr

from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.spin_mecp import MECPConfig, SurfacePair, optimize_mecp, seam_curvature
from pincer_catmech.quantum.spin_provider import validate_spin_system


REPO = Path(__file__).resolve().parents[1]
PROVIDER_SOURCE = REPO / "src/pincer_catmech/quantum/spin_provider.py"
SETTING_NAMES = ("method", "basis", "scf_max_iterations", "grid_radial", "grid_spherical",
                 "scf_algorithm", "e_convergence", "d_convergence")


def resolve_evidence_path(path, *, relative_to=None):
    """Resolve read access only; never rewrite recorded paths or signatures.

    Existing literal paths retain priority. A missing historical Windows/POSIX
    path may relocate only its exact ``data/`` suffix into this checkout's
    data directory. Ambiguous matches and escaping/symlinked destinations fail
    closed. Engine binaries and other missing external paths are not relocated.
    """
    original = Path(path)
    if original.exists():
        return original
    relative = original
    if (relative_to is not None and not original.is_absolute()
            and not PureWindowsPath(str(path)).is_absolute()):
        relative = Path(relative_to)/original
        if relative.exists():
            return relative
    parts = str(path).replace("\\", "/").split("/")
    data_root = REPO.resolve()/"data"
    candidates = []
    for index, part in enumerate(parts):
        if part != "data" or index == len(parts)-1:
            continue
        suffix = parts[index+1:]
        _require(not any(piece in ("", ".", "..") for piece in suffix), "Unsafe historical data path")
        candidate = data_root.joinpath(*suffix).resolve()
        _require(candidate.is_relative_to(data_root), "Historical data path escapes repository data directory")
        if candidate.is_file() and candidate not in candidates:
            candidates.append(candidate)
    _require(len(candidates) <= 1, "Ambiguous historical data path relocation")
    return candidates[0] if candidates else relative


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def geometry_digest(symbols, positions):
    return hashlib.sha256(np.asarray(positions, dtype="<f8").tobytes()
                          + " ".join(symbols).encode()).hexdigest()


def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _native_protocol(response, symbols, positions, charge, protocol, multiplicity):
    _require(response["geometry_sha256"] == geometry_digest(symbols, positions), "native_geometry_hash_mismatch")
    _require(response["symbols"] == symbols and response["charge"] == charge, "native_identity_mismatch")
    _require(np.shape(response["positions_bohr"]) == positions.shape and np.allclose(
        response["positions_bohr"], positions, atol=1e-10, rtol=0), "native_geometry_mismatch")
    _require(response["multiplicities"] == [multiplicity], "native_multiplicity_mismatch")
    for key, expected in (("coordinate_units", "bohr"), ("energy_units", "hartree"),
                          ("gradient_units", "hartree/bohr")):
        _require(response.get(key) == expected, f"native_units_mismatch:{key}")
    for key in ("method", "basis", "grid_radial", "grid_spherical", "scf_algorithm", "scf_max_iterations"):
        _require(response[key] == protocol[key], f"native_protocol_mismatch:{key}")
    for native_key, key in (("scf_e_convergence", "e_convergence"), ("scf_d_convergence", "d_convergence")):
        _require(response[native_key] == protocol[key], f"native_protocol_mismatch:{key}")
    if protocol["scf_algorithm"] == "soscf":
        _require(response["soscf_start_convergence"] == protocol["soscf_start_convergence"], "native_SOSCF_threshold_mismatch")


def _verified_state(row, python_executable, identity_checker):
    """Recheck original evidence, numeric labels and current structural identity."""
    _require(row.get("status") == "converged", "state_not_converged_or_missing")
    _require(row.get("scf_converged") is True, "SCF_quality_not_passed")
    _require(row.get("intended_catalyst_identity_pass") is True, "identity_quality_not_passed")
    _require(row.get("spin_contamination_flag") is False, "spin_quality_not_passed")
    _require(row.get("gradient_kind") == "analytic", "analytic_gradient_required")
    symbols, positions, multiplicity = row["symbols"], np.asarray(row["positions_bohr"], dtype=float), row["multiplicity"]
    validate_spin_system(symbols, positions, row["charge"], [multiplicity], require_pair=False)
    expected = (multiplicity*multiplicity-1)/4
    _require(np.ndim(row["spin_squared"]) == 0 and np.isfinite(row["spin_squared"])
             and abs(row["spin_squared"]-expected) <= .1, "S2_eligibility_screen_failed")
    _require(row["expected_spin_squared"] == expected, "S2_reference_mismatch")
    _require(np.isfinite(row["energy_hartree"]), "nonfinite_energy")
    gradient = np.asarray(row["gradient_hartree_bohr"], dtype=float)
    _require(gradient.shape == positions.shape and np.all(np.isfinite(gradient)), "invalid_gradient")
    _require(geometry_digest(symbols, positions) == row["geometry_sha256"], "geometry_hash_mismatch")
    protocol = row["protocol"]
    _require(isinstance(protocol, dict) and isinstance(protocol.get("method"), str), "malformed_seed_protocol")
    _require(protocol.get("provider") == "psi4" and protocol.get("scf_type") == "df"
             and protocol.get("solvent") is None, "unsupported_seed_protocol")
    _require(json_digest(protocol) == row["protocol_sha256"], "protocol_hash_mismatch")
    paired = {key: value for key, value in protocol.items() if key != "reference"}
    _require(json_digest(paired) == row["paired_protocol_sha256"], "paired_protocol_hash_mismatch")
    _require(protocol["provider_source_sha256"] == digest(PROVIDER_SOURCE), "provider_source_changed")
    _require(protocol["executable_sha256"] == digest(python_executable), "Psi4_executable_changed")
    fingerprint = protocol.get("engine_fingerprint")
    _require(isinstance(fingerprint, dict) and bool(fingerprint), "missing_engine_fingerprint")
    prefix = Path(python_executable).resolve().parent
    for relative, checksum in fingerprint.items():
        target = (prefix/relative).resolve()
        _require(target.is_relative_to(prefix), "engine_fingerprint_path_outside_environment")
        _require(digest(target) == checksum, "Psi4_engine_fingerprint_changed")
    for name in SETTING_NAMES:
        _require(name in protocol, f"missing_protocol_setting:{name}")
    _require(protocol["scf_algorithm"] in ("diis", "soscf"), "unsupported_SCF_algorithm")
    if protocol["scf_algorithm"] == "soscf":
        _require(protocol.get("soscf_max_iter") == 5, "unsupported_SOSCF_protocol")
        _require(np.isfinite(protocol["soscf_start_convergence"])
                 and protocol["soscf_start_convergence"] > 0, "invalid_SOSCF_threshold")
    hf = protocol["method"].lower() in ("hf", "scf")
    reference = ("rhf" if multiplicity == 1 else "uhf") if hf else ("rks" if multiplicity == 1 else "uks")
    _require(protocol.get("reference") == reference, "state_reference_mismatch")
    signature = dict(symbols=symbols, positions_bohr=positions.tolist(), charge=row["charge"],
                     multiplicity=multiplicity, protocol=protocol)
    _require(json_digest(signature) == row["signature_sha256"], "state_signature_mismatch")
    for path_key, hash_key in (("source_xyz", "source_xyz_sha256"),
                               ("source_metadata_path", "source_metadata_sha256")):
        _require(digest(resolve_evidence_path(row[path_key])) == row[hash_key], f"source_changed:{path_key}")
    atoms = read(resolve_evidence_path(row["source_xyz"]))
    _require(atoms.get_chemical_symbols() == symbols and np.allclose(
        atoms.positions/Bohr, positions, atol=1e-10, rtol=0), "source_XYZ_geometry_mismatch")
    metadata = json.loads(resolve_evidence_path(row["source_metadata_path"]).read_text(encoding="utf-8"))
    _require(metadata["charge"] == row["charge"], "source_charge_mismatch")
    _require(identity_checker(atoms, metadata).get("retained") is True, "current_identity_screen_failed")
    evidence = row.get("native_evidence", [])
    _require(bool(evidence), "missing_native_evidence")
    responses = []
    for item in evidence:
        _require(digest(resolve_evidence_path(item["path"])) == item["sha256"], "native_evidence_hash_mismatch")
        if resolve_evidence_path(item["path"]).name == "response.json":
            responses.append(item["path"])
    _require(len(responses) == 1, "one_native_state_response_required")
    native = json.loads(resolve_evidence_path(responses[0]).read_text(encoding="utf-8"))
    _native_protocol(native, symbols, positions, row["charge"], protocol, multiplicity)
    _require(native["geometry_sha256"] == row["geometry_sha256"] and len(native["states"]) == 1,
             "native_geometry_or_state_count_mismatch")
    state = native["states"][0]
    for key in ("multiplicity", "energy_hartree", "gradient_hartree_bohr", "spin_squared",
                "expected_spin_squared", "scf_converged", "gradient_kind"):
        _require(state[key] == row[key], f"native_numeric_mismatch:{key}")
    return metadata


def inspect_matrix(matrix_path, python_executable, *, catalyst=None, identity_checker=catalyst_identity_screen):
    """Read-only eligibility gate. Incomplete matrices never launch MECP work."""
    matrix_path = resolve_evidence_path(matrix_path).resolve()
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    _require(isinstance(matrix, dict), "matrix_summary_must_be_an_object")
    if "active_matrix_summary" in matrix:
        child = resolve_evidence_path(matrix["active_matrix_summary"], relative_to=matrix_path.parent)
        _require(digest(child) == matrix["active_matrix_summary_sha256"], "active_matrix_summary_hash_mismatch")
        matrix_path = child.resolve()
        matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
        _require(isinstance(matrix, dict), "matrix_summary_must_be_an_object")
    rows = matrix.get("state_records", [])
    ready = dict(matrix_summary=str(matrix_path), matrix_sha256=digest(matrix_path),
                 matrix_run_complete=matrix.get("run_complete") is True,
                 planned_target_states=matrix.get("planned_target_states"), state_records=len(rows),
                 eligible_pairs=[], rejected_states=[], selected_seed=None)
    complete = (matrix.get("run_complete") is True and isinstance(rows, list)
                and isinstance(matrix.get("planned_target_states"), int)
                and not isinstance(matrix.get("planned_target_states"), bool)
                and matrix["planned_target_states"] > 0
                and len(rows) == matrix["planned_target_states"] == matrix.get("target_states"))
    if not complete:
        return dict(ready, status="matrix_incomplete", reason="Wait for every planned state or explicit missing-state record")
    _require(all(isinstance(row, dict) for row in rows), "malformed_matrix_state_record")
    identities = [(row.get("catalyst_id"), row.get("multiplicity")) for row in rows]
    _require(len(identities) == len(set(identities)), "duplicate_matrix_state_identity")
    valid, source_metadata = [], {}
    for row in rows:
        if catalyst is not None and row.get("catalyst_id") != catalyst:
            continue
        try:
            _require(row.get("chemical_state") == "active", "active_catalyst_seed_required")
            _require(isinstance(row.get("catalyst_id"), str)
                     and row["catalyst_id"].startswith(("Fe_", "Co_", "Mn_")), "outside_target_spin_scope")
            metadata = _verified_state(row, Path(python_executable), identity_checker)
            valid.append(row)
            source_metadata[row["signature_sha256"]] = metadata
        except (KeyError, TypeError, ValueError, OSError) as exc:
            ready["rejected_states"].append(dict(catalyst_id=row.get("catalyst_id"),
                multiplicity=row.get("multiplicity"), reason=str(exc)))
    seeds = []
    for a, b in itertools.combinations(valid, 2):
        keys = ("catalyst_id", "geometry_sha256", "paired_protocol_sha256", "charge", "source_metadata_sha256")
        if any(a[key] != b[key] for key in keys) or a["multiplicity"] == b["multiplicity"]:
            continue
        a, b = sorted((a, b), key=lambda state: state["multiplicity"])
        pair = SurfacePair(a["energy_hartree"], b["energy_hartree"],
                           np.asarray(a["gradient_hartree_bohr"]), np.asarray(b["gradient_hartree_bohr"]),
                           f"multiplicity_{a['multiplicity']}", f"multiplicity_{b['multiplicity']}", True,
                           dict(states=[a, b], intended_catalyst_identity_pass=True))
        pair.validate(np.shape(a["positions_bohr"]), positions_bohr=np.asarray(a["positions_bohr"]))
        seed = dict(catalyst_id=a["catalyst_id"], multiplicities=[a["multiplicity"], b["multiplicity"]],
                    charge=a["charge"], symbols=a["symbols"], positions_bohr=a["positions_bohr"],
                    geometry_sha256=a["geometry_sha256"], paired_protocol_sha256=a["paired_protocol_sha256"],
                    source_state_signatures=[a["signature_sha256"], b["signature_sha256"]],
                    source_metadata=source_metadata[a["signature_sha256"]], protocol=a["protocol"],
                    vertical_gap_hartree=pair.gap)
        seeds.append(seed)
    seeds.sort(key=lambda seed: (abs(seed["vertical_gap_hartree"]), seed["catalyst_id"], seed["multiplicities"]))
    ready["eligible_pairs"] = [{key: seed[key] for key in ("catalyst_id", "multiplicities",
        "vertical_gap_hartree", "geometry_sha256", "paired_protocol_sha256")} for seed in seeds]
    if not seeds:
        return dict(ready, status="no_eligible_pair", reason="No two currently verified states pass all seed gates at matched geometry/protocol")
    return dict(ready, status="eligible_pair_ready", selected_seed=seeds[0],
                selection_rule="smallest absolute vertical gap among eligible matched pairs; no crossing inference from gap sign")


def _kill_tree(process, measured):
    try:
        descendants = measured.children(recursive=True)
    except psutil.NoSuchProcess:
        descendants = []
    for child in reversed(descendants):
        try:
            child.kill()
        except psutil.NoSuchProcess:
            pass
    if process.poll() is None:
        process.kill()
    process.wait()


def run_worker(command, *, cwd, env, output, timeout_seconds, rss_limit_mib=600):
    """Monitor the owned worker tree; reject late/over-memory results."""
    start, peak, samples, status, code = time.monotonic(), 0, 0, "completed", None
    with Path(output).open("w", encoding="utf-8") as log:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        try:
            measured = psutil.Process(process.pid)
        except psutil.NoSuchProcess:
            code = process.wait()
            elapsed = time.monotonic()-start
            return dict(status="timeout" if elapsed > timeout_seconds else ("memory_monitor_unavailable" if code == 0 else "backend_failure"),
                        returncode=code, wall_seconds=elapsed, peak_rss_mib=None, rss_limit_mib=rss_limit_mib,
                        timeout_seconds=timeout_seconds, argv=command, memory_sampling="worker exited before first sample")
        except psutil.Error:
            if process.poll() is None:
                process.kill()
            process.wait()
            raise
        try:
            while process.poll() is None:
                try:
                    descendants = measured.children(recursive=True)
                    rss = measured.memory_info().rss + sum(child.memory_info().rss for child in descendants)
                    peak = max(peak, rss)
                    samples += 1
                except psutil.NoSuchProcess:
                    if process.poll() is not None:
                        break
                if peak > rss_limit_mib*1024**2:
                    status = "memory_limit"
                    break
                if time.monotonic()-start >= timeout_seconds:
                    status = "timeout"
                    break
                time.sleep(.1)
        finally:
            _kill_tree(process, measured)
            code = process.returncode
    elapsed = time.monotonic()-start
    if status == "completed" and elapsed > timeout_seconds:
        status = "timeout"
    if status == "completed" and code:
        status = "backend_failure"
    if status == "completed" and samples == 0:
        status = "memory_monitor_unavailable"
    return dict(status=status, returncode=code, wall_seconds=elapsed, peak_rss_mib=peak/1024**2,
                rss_samples=samples, rss_limit_mib=rss_limit_mib, timeout_seconds=timeout_seconds, argv=command)


class BoundedPairProvider:
    """Two serial 120-second state workers; one shared pair/total-wall budget."""
    def __init__(self, seed, python_executable, workdir, *, max_pair_evaluations=8,
                 wall_seconds=1920, identity_checker=catalyst_identity_screen):
        _require(isinstance(max_pair_evaluations, int) and not isinstance(max_pair_evaluations, bool)
                 and max_pair_evaluations > 0, "Positive pair evaluation budget required")
        _require(np.isfinite(wall_seconds) and wall_seconds > 0, "Positive finite wall budget required")
        self.seed, self.python = seed, Path(python_executable).resolve()
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.limit, self.wall_seconds, self.started = max_pair_evaluations, wall_seconds, time.monotonic()
        self.evaluations, self.state_attempts = 0, 0
        self.identity_checker = identity_checker
        self.records = []

    def __call__(self, positions_bohr):
        remaining = self.wall_seconds-(time.monotonic()-self.started)
        if self.evaluations >= self.limit or remaining <= 0:
            raise TimeoutError("Shared MECP/curvature pair or wall budget exhausted")
        seed = self.seed
        x = validate_spin_system(seed["symbols"], positions_bohr, seed["charge"], seed["multiplicities"])
        identity = self.identity_checker(Atoms(seed["symbols"], positions=x*Bohr), seed["source_metadata"])
        _require(identity.get("retained") is True, "Trial geometry failed catalyst identity screen")
        _require(digest(PROVIDER_SOURCE) == seed["protocol"]["provider_source_sha256"], "Provider changed during MECP")
        _require(digest(self.python) == seed["protocol"]["executable_sha256"], "Psi4 executable changed during MECP")
        for relative, checksum in seed["protocol"]["engine_fingerprint"].items():
            _require(digest(self.python.parent/relative) == checksum, "Psi4 engine changed during MECP")
        self.evaluations += 1
        pair_dir = self.workdir/f"pair_{self.evaluations:04d}_{uuid.uuid4().hex[:8]}"
        pair_dir.mkdir()
        states = []
        protocol = seed["protocol"]
        settings = {key: protocol[key] for key in SETTING_NAMES}
        settings["soscf_start_convergence"] = protocol.get("soscf_start_convergence") or 1e-4
        env = dict(os.environ)
        env.update(OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
                   PYTHONUTF8="1", PYTHONIOENCODING="utf-8", PYTHONPATH=str(REPO/"src"))
        prefix = self.python.parent
        env["PATH"] = str(prefix/"Library/bin")+os.pathsep+str(prefix)+os.pathsep+env.get("PATH", "")
        for multiplicity in seed["multiplicities"]:
            remaining = self.wall_seconds-(time.monotonic()-self.started)
            if remaining <= 0:
                raise TimeoutError("Shared MECP wall budget exhausted before next state")
            folder = pair_dir/f"multiplicity_{multiplicity}"
            folder.mkdir()
            request = dict(symbols=seed["symbols"], positions_bohr=x.tolist(), charge=seed["charge"],
                           multiplicities=[multiplicity], single_state=True, **settings,
                           memory_mib=500, threads=2, workdir=str(folder/"native"))
            dump(folder/"request.json", request)
            command = [str(self.python), "-m", "pincer_catmech.quantum.spin_provider", "--request",
                       str(folder/"request.json"), "--result", str(folder/"response.json")]
            self.state_attempts += 1
            try:
                launch = run_worker(command, cwd=folder, env=env, output=folder/"worker.out",
                                    timeout_seconds=min(120., remaining), rss_limit_mib=600)
            except (OSError, psutil.Error) as exc:
                launch = dict(status="launch_or_monitor_failure", failure=str(exc), argv=command)
            launch.update(multiplicity=multiplicity, pair_evaluation=self.evaluations, configured_memory_mib=500)
            dump(folder/"launch.json", launch)
            self.records.append(dict(path=str(folder/"launch.json"), **launch))
            if launch["status"] != "completed":
                failure = f"MECP state {multiplicity}: {launch['status']}; inspect {folder}"
                if launch["status"] == "timeout":
                    raise TimeoutError(failure)
                raise RuntimeError(failure)
            response = json.loads((folder/"response.json").read_text(encoding="utf-8"))
            _native_protocol(response, seed["symbols"], x, seed["charge"], protocol, multiplicity)
            _require(len(response["states"]) == 1 and response["states"][0]["multiplicity"] == multiplicity,
                     "Worker returned wrong electronic state")
            state = response["states"][0]
            _require(state.get("scf_converged") is True and state.get("gradient_kind") == "analytic", "Worker state quality failed")
            expected = (multiplicity*multiplicity-1)/4
            _require(state["expected_spin_squared"] == expected and np.isfinite(state["spin_squared"])
                     and abs(state["spin_squared"]-expected) <= .1, "Worker S2 eligibility screen failed")
            states.append(state)
        if time.monotonic()-self.started > self.wall_seconds:
            raise TimeoutError("MECP wall budget exceeded by completed pair")
        a, b = states
        metadata = dict(states=states, positions_bohr=x.tolist(), symbols=seed["symbols"], charge=seed["charge"],
                        multiplicities=seed["multiplicities"], **settings, coordinate_units="bohr",
                        energy_units="hartree", gradient_units="hartree/bohr", intended_catalyst_identity_pass=True,
                        identity_validation=identity, paired_protocol_sha256=seed["paired_protocol_sha256"])
        pair = SurfacePair(a["energy_hartree"], b["energy_hartree"], np.asarray(a["gradient_hartree_bohr"]),
                           np.asarray(b["gradient_hartree_bohr"]), f"multiplicity_{a['multiplicity']}",
                           f"multiplicity_{b['multiplicity']}", True, metadata)
        pair.validate(x.shape, positions_bohr=x)
        dump(pair_dir/"pair.json", dict(energy1_hartree=pair.energy1_hartree, energy2_hartree=pair.energy2_hartree,
                                      gap_hartree=pair.gap, metadata=metadata))
        return pair


def run(matrix_path, python_executable, output_dir, *, max_pair_evaluations=8, wall_seconds=1920,
        readiness_only=False, catalyst=None, provider_factory=BoundedPairProvider,
        identity_checker=catalyst_identity_screen):
    """Write an immutable attempt plus a latest-summary pointer; never alter input."""
    _require(isinstance(max_pair_evaluations, int) and not isinstance(max_pair_evaluations, bool)
             and max_pair_evaluations > 0, "Positive pair evaluation budget required")
    _require(np.isfinite(wall_seconds) and wall_seconds > 0, "Positive finite wall budget required")
    output_dir = Path(output_dir).resolve()
    attempt = output_dir/f"attempt_{uuid.uuid4().hex}"
    attempt.mkdir(parents=True, exist_ok=False)
    receipt = dict(started_utc=datetime.now(timezone.utc).isoformat(), quantum_pair_evaluations=0,
                   quantum_state_attempts=0, max_pair_evaluations=max_pair_evaluations,
                   wall_budget_seconds=wall_seconds, per_state_timeout_seconds=120,
                   configured_memory_mib=500, observed_rss_limit_mib=600, threads=2,
                   first_order_crossing=False, minimum_verified=False, validated_catalyst_label=False,
                   provider_source_sha256=digest(PROVIDER_SOURCE), runner_source_sha256=digest(__file__),
                   optimizer_source_sha256=digest(REPO/"src/pincer_catmech/quantum/spin_mecp.py"),
                   evidence="prospective bounded pincer crossing search; no SOC, rate, experiment or global minimum claim")
    try:
        readiness = inspect_matrix(matrix_path, python_executable, catalyst=catalyst, identity_checker=identity_checker)
        receipt.update(readiness=readiness, status=readiness["status"])
        if readiness["status"] == "eligible_pair_ready" and not readiness_only:
            seed = readiness["selected_seed"]
            provider = provider_factory(seed, python_executable, attempt/"native", max_pair_evaluations=max_pair_evaluations,
                                        wall_seconds=wall_seconds, identity_checker=identity_checker)
            cfg = MECPConfig(max_evaluations=max_pair_evaluations, wall_seconds=wall_seconds)
            result = optimize_mecp(provider, seed["positions_bohr"], cfg)
            receipt.update(status=result.status, optimizer=result.to_dict(), first_order_crossing=result.converged)
            if result.converged:
                try:
                    curvature = seam_curvature(provider, result.positions_bohr)
                    receipt.update(curvature=curvature, minimum_verified=curvature["minimum_verified"],
                                   status="local_numerical_minimum" if curvature["minimum_verified"] else "first_order_crossing_curvature_not_accepted")
                except (TimeoutError, ValueError, RuntimeError, OSError) as exc:
                    receipt.update(status="first_order_crossing_curvature_incomplete", curvature_failure=str(exc))
            receipt.update(quantum_pair_evaluations=provider.evaluations, quantum_state_attempts=provider.state_attempts,
                           worker_records=provider.records)
    except (KeyError, TypeError, ValueError, RuntimeError, OSError) as exc:
        receipt.update(status="invalid_matrix_or_orchestration_failure", failure=str(exc))
    evidence = [dict(path=str(path), sha256=digest(path)) for path in sorted(attempt.rglob("*")) if path.is_file()]
    receipt.update(completed_utc=datetime.now(timezone.utc).isoformat(), native_evidence=evidence)
    dump(attempt/"result.json", receipt)
    dump(output_dir/"summary.json", dict(**receipt, attempt_result=str(attempt/"result.json"),
                                        attempt_result_sha256=digest(attempt/"result.json")))
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-summary", type=Path, required=True)
    parser.add_argument("--python", "--psi4-python", dest="python", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("data/phase3/spin/pincer_mecp_from_matrix_001"))
    parser.add_argument("--max-pair-evaluations", type=int, default=8)
    parser.add_argument("--wall-seconds", type=float, default=1920)
    parser.add_argument("--readiness-only", action="store_true")
    parser.add_argument("--catalyst")
    args = parser.parse_args(argv)
    receipt = run(args.matrix_summary, args.python, args.output_dir, max_pair_evaluations=args.max_pair_evaluations,
                  wall_seconds=args.wall_seconds, readiness_only=args.readiness_only, catalyst=args.catalyst)
    print(json.dumps({key: receipt[key] for key in ("status", "quantum_pair_evaluations", "quantum_state_attempts",
                                                  "first_order_crossing", "minimum_verified")}))
    return 2 if receipt["status"] in ("matrix_incomplete", "invalid_matrix_or_orchestration_failure") else 0


if __name__ == "__main__":
    sys.exit(main())
