"""Prospective orchestration tests: every electronic-structure worker is mocked."""
from copy import deepcopy
import importlib.util
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write
from ase.units import Bohr


@pytest.fixture
def cli():
    spec = importlib.util.spec_from_file_location("phase3_mecp_cli", Path(__file__).resolve().parents[1]/"scripts/run_phase3_mecp.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def identity_pass(*args):
    return {"retained": True, "evidence": "mock identity for orchestration tests only"}


def native_response(cli, seed, multiplicity, positions, protocol):
    x = np.asarray(positions)
    delta = x-np.asarray(seed["positions_bohr"])
    sign = 1 if multiplicity == 1 else -1
    normal = np.zeros_like(delta)
    normal[0, 2], normal[1, 2] = 1., -1.
    gradient = delta+sign*normal
    expected = (multiplicity*multiplicity-1)/4
    state = dict(multiplicity=multiplicity, energy_hartree=float(np.sum(delta*delta)/2+sign*np.sum(delta*normal)+.5),
                 gradient_hartree_bohr=gradient.tolist(), spin_squared=expected, expected_spin_squared=expected,
                 scf_converged=True, gradient_kind="analytic", reference="rks" if multiplicity == 1 else "uks")
    return dict(symbols=seed["symbols"], positions_bohr=x.tolist(), charge=0, multiplicities=[multiplicity],
                geometry_sha256=cli.geometry_digest(seed["symbols"], x), states=[state], coordinate_units="bohr",
                energy_units="hartree", gradient_units="hartree/bohr", psi4_version="mock",
                **{key:protocol[key] for key in ("method", "basis", "grid_radial", "grid_spherical",
                    "scf_algorithm", "scf_max_iterations", "soscf_start_convergence")},
                scf_e_convergence=protocol["e_convergence"], scf_d_convergence=protocol["d_convergence"])


@pytest.fixture
def matrix(tmp_path, cli):
    engine = tmp_path/"engine"
    engine.mkdir()
    python = engine/"python.exe"
    python.write_bytes(b"never execute this test fixture")
    (engine/"engine.bin").write_bytes(b"mock engine fingerprint")
    xyz = tmp_path/"best_found.xyz"
    write(xyz, Atoms("FeHH", positions=[[0.,0.,0.], [0.,0.,2.], [0.,2.,0.]]))
    metadata = tmp_path/"metadata.json"
    cli.dump(metadata, {"charge":0, "state":"active"})
    atoms = read(xyz)
    seed = dict(symbols=atoms.get_chemical_symbols(), positions_bohr=(atoms.positions/Bohr).tolist())
    protocol = dict(provider="psi4", method="pbe", basis="sto-3g", solvent=None, scf_type="df",
                    e_convergence=1e-8, d_convergence=1e-6, scf_max_iterations=120, grid_radial=35,
                    grid_spherical=110, scf_algorithm="soscf", soscf_start_convergence=1e-4, soscf_max_iter=5,
                    engine_fingerprint={"engine.bin":cli.digest(engine/"engine.bin")},
                    executable_sha256=cli.digest(python), provider_source_sha256=cli.digest(cli.PROVIDER_SOURCE))
    rows = []
    for multiplicity in (1,3):
        state_protocol = dict(protocol, reference="rks" if multiplicity == 1 else "uks")
        native = native_response(cli, seed, multiplicity, seed["positions_bohr"], state_protocol)
        folder = tmp_path/f"state_{multiplicity}"
        folder.mkdir()
        response = folder/"response.json"
        cli.dump(response, native)
        signature = dict(**seed, charge=0, multiplicity=multiplicity, protocol=state_protocol)
        rows.append(dict(catalyst_id="Fe_mock_fixture", chemical_state="active", charge=0,
                         **seed, **native["states"][0], status="converged", spin_contamination_flag=False,
                         intended_catalyst_identity_pass=True, geometry_sha256=native["geometry_sha256"],
                         protocol=state_protocol, protocol_sha256=cli.json_digest(state_protocol),
                         paired_protocol_sha256=cli.json_digest(protocol), signature_sha256=cli.json_digest(signature),
                         source_xyz=str(xyz), source_xyz_sha256=cli.digest(xyz), source_metadata_path=str(metadata),
                         source_metadata_sha256=cli.digest(metadata), native_evidence=[dict(path=str(response),sha256=cli.digest(response))]))
    payload = dict(run_complete=True, planned_target_states=2, target_states=2, state_records=rows)
    path = tmp_path/"matrix.json"
    cli.dump(path, payload)
    return SimpleNamespace(path=path, python=python, payload=payload, seed=seed, protocol=protocol, root=tmp_path)


def test_completed_matrix_can_propose_seed_without_claiming_crossing(cli, matrix):
    ready = cli.inspect_matrix(matrix.path, matrix.python, identity_checker=identity_pass)
    assert ready["status"] == "eligible_pair_ready"
    assert ready["selected_seed"]["multiplicities"] == [1,3]
    assert len(ready["eligible_pairs"]) == 1
    assert "no crossing inference" in ready["selection_rule"]
    receipt = cli.run(matrix.path, matrix.python, matrix.root/"out", readiness_only=True,
                      identity_checker=identity_pass)
    assert receipt["status"] == "eligible_pair_ready"
    assert receipt["quantum_state_attempts"] == 0
    assert not receipt["first_order_crossing"] and not receipt["minimum_verified"]


@pytest.mark.parametrize("edit", [
    {"spin_contamination_flag":True}, {"spin_squared":2.6614}, {"scf_converged":False},
    {"intended_catalyst_identity_pass":False}, {"status":"timeout"},
    {"energy_hartree":99.}, {"expected_spin_squared":3.},
])
def test_failed_or_tampered_state_produces_zero_call_no_pair_receipt(cli, matrix, edit):
    matrix.payload["state_records"][1].update(edit)
    cli.dump(matrix.path, matrix.payload)
    def forbidden(*args, **kwargs):
        pytest.fail("An ineligible matrix must not create a quantum provider")
    result = cli.run(matrix.path, matrix.python, matrix.root/"out", provider_factory=forbidden,
                     identity_checker=identity_pass)
    assert result["status"] == "no_eligible_pair"
    assert result["quantum_pair_evaluations"] == result["quantum_state_attempts"] == 0
    assert not result["minimum_verified"] and result["readiness"]["rejected_states"]
    assert (matrix.root/"out/summary.json").is_file()


def test_incomplete_matrix_waits_even_if_recorded_rows_are_eligible(cli, matrix):
    matrix.payload.update(run_complete=False, planned_target_states=3)
    cli.dump(matrix.path, matrix.payload)
    ready = cli.inspect_matrix(matrix.path, matrix.python, identity_checker=identity_pass)
    assert ready["status"] == "matrix_incomplete" and ready["selected_seed"] is None


@pytest.mark.parametrize("target", ["native", "engine", "geometry"])
def test_original_evidence_must_still_match_its_hash(cli, matrix, target):
    row = matrix.payload["state_records"][1]
    path = {"native":Path(row["native_evidence"][0]["path"]),
            "engine":matrix.python.parent/"engine.bin", "geometry":Path(row["source_xyz"])}[target]
    path.write_bytes(path.read_bytes()+b"changed")
    ready = cli.inspect_matrix(matrix.path, matrix.python, identity_checker=identity_pass)
    assert ready["status"] == "no_eligible_pair"


def test_current_identity_screen_overrides_old_pass_flag(cli, matrix):
    ready = cli.inspect_matrix(matrix.path, matrix.python, identity_checker=lambda *args:{"retained":False})
    assert ready["status"] == "no_eligible_pair"
    assert all("current_identity" in row["reason"] for row in ready["rejected_states"])


def test_native_protocol_cannot_be_relabelled_by_changing_row_hashes(cli, matrix):
    row = matrix.payload["state_records"][1]
    row["protocol"]["basis"] = "def2-svp"
    row["protocol_sha256"] = cli.json_digest(row["protocol"])
    row["paired_protocol_sha256"] = cli.json_digest({k:v for k,v in row["protocol"].items() if k!="reference"})
    row["signature_sha256"] = cli.json_digest(dict(symbols=row["symbols"], positions_bohr=row["positions_bohr"],
        charge=row["charge"], multiplicity=row["multiplicity"], protocol=row["protocol"]))
    cli.dump(matrix.path, matrix.payload)
    ready = cli.inspect_matrix(matrix.path, matrix.python, identity_checker=identity_pass)
    assert ready["status"] == "no_eligible_pair"
    assert "native_protocol_mismatch:basis" in ready["rejected_states"][0]["reason"]


def test_aggregate_pointer_requires_matching_child_hash(cli, matrix):
    aggregate = matrix.root/"aggregate.json"
    cli.dump(aggregate, dict(active_matrix_summary=str(matrix.path), active_matrix_summary_sha256=cli.digest(matrix.path)))
    assert cli.inspect_matrix(aggregate, matrix.python, identity_checker=identity_pass)["status"] == "eligible_pair_ready"
    matrix.path.write_text(matrix.path.read_text()+" ")
    with pytest.raises(ValueError, match="active_matrix_summary_hash_mismatch"):
        cli.inspect_matrix(aggregate, matrix.python, identity_checker=identity_pass)


def install_mock_worker(monkeypatch, cli, matrix, *, failure=None, tangent_offset=False):
    requests = []
    def worker(command, **kwargs):
        request = json.loads(Path(command[command.index("--request")+1]).read_text())
        requests.append(request)
        assert request["memory_mib"] == 500 and request["threads"] == 2
        assert 0 < kwargs["timeout_seconds"] <= 120 and kwargs["rss_limit_mib"] == 600
        response = native_response(cli, matrix.seed, request["multiplicities"][0], request["positions_bohr"], matrix.protocol)
        if failure == "spin":
            response["states"][0]["spin_squared"] += .8
        elif failure == "protocol":
            response["basis"] = "unexpected"
        elif failure == "units":
            response["coordinate_units"] = "angstrom"
        if tangent_offset:
            response["states"][0]["gradient_hartree_bohr"][1][1] += .2
            response["states"][0]["energy_hartree"] += .2*(request["positions_bohr"][1][1]-matrix.seed["positions_bohr"][1][1])
        cli.dump(command[command.index("--result")+1], response)
        return dict(status="completed", returncode=0, wall_seconds=.001, peak_rss_mib=10.)
    monkeypatch.setattr(cli, "run_worker", worker)
    return requests


def test_search_and_curvature_share_pair_budget_and_do_not_fake_minimum(cli, matrix, monkeypatch):
    requests = install_mock_worker(monkeypatch, cli, matrix)
    receipt = cli.run(matrix.path, matrix.python, matrix.root/"out", max_pair_evaluations=1,
                      identity_checker=identity_pass)
    assert receipt["first_order_crossing"]
    assert receipt["status"] == "first_order_crossing_curvature_incomplete"
    assert not receipt["minimum_verified"]
    assert receipt["quantum_pair_evaluations"] == 1 and len(requests) == receipt["quantum_state_attempts"] == 2
    summary = json.loads((matrix.root/"out/summary.json").read_text())
    assert cli.digest(summary["attempt_result"]) == summary["attempt_result_sha256"]


def test_analytic_mock_candidate_requires_real_curvature_calls_in_budget(cli, matrix, monkeypatch):
    requests = install_mock_worker(monkeypatch, cli, matrix)
    receipt = cli.run(matrix.path, matrix.python, matrix.root/"out", max_pair_evaluations=8,
                      identity_checker=identity_pass)
    assert receipt["status"] == "local_numerical_minimum", receipt
    assert receipt["minimum_verified"] and not receipt["validated_catalyst_label"]
    assert receipt["quantum_pair_evaluations"] == 6
    assert receipt["quantum_state_attempts"] == len(requests) == 12
    assert receipt["curvature"]["hessian_symmetry_pass"]


@pytest.mark.parametrize("failure", ["spin", "protocol", "units"])
def test_new_worker_quality_failure_stops_before_second_state(cli, matrix, monkeypatch, failure):
    requests = install_mock_worker(monkeypatch, cli, matrix, failure=failure)
    receipt = cli.run(matrix.path, matrix.python, matrix.root/"out", identity_checker=identity_pass)
    assert receipt["status"] == "invalid_surfaces_or_backend_failure"
    assert not receipt["first_order_crossing"] and not receipt["minimum_verified"]
    assert len(requests) == receipt["quantum_state_attempts"] == 1


def test_nonstationary_optimizer_exhausts_exact_shared_budget(cli, matrix, monkeypatch):
    requests = install_mock_worker(monkeypatch, cli, matrix, tangent_offset=True)
    receipt = cli.run(matrix.path, matrix.python, matrix.root/"out", max_pair_evaluations=1,
                      identity_checker=identity_pass)
    assert receipt["status"] == "budget_exhausted"
    assert not receipt["first_order_crossing"] and not receipt["minimum_verified"]
    assert receipt["quantum_pair_evaluations"] == 1 and len(requests) == 2


@pytest.mark.parametrize("limit", ["memory", "timeout"])
def test_worker_monitor_terminates_owned_process_on_resource_limit(cli, matrix, monkeypatch, limit):
    class Process:
        pid = 123
        returncode = None
        killed = False
        def poll(self):
            return self.returncode
        def kill(self):
            self.killed = True
            self.returncode = -9
        def wait(self):
            return self.returncode
    process = Process()
    measured = SimpleNamespace(children=lambda **kwargs:[], memory_info=lambda:SimpleNamespace(rss=(601 if limit=="memory" else 10)*1024**2))
    monkeypatch.setattr(cli.subprocess, "Popen", lambda *args, **kwargs:process)
    monkeypatch.setattr(cli.psutil, "Process", lambda pid:measured)
    times = iter([0., 121., 121.])
    monkeypatch.setattr(cli.time, "monotonic", lambda:next(times, 121.))
    result = cli.run_worker(["never-executed"], cwd=matrix.root, env={}, output=matrix.root/"worker.out", timeout_seconds=120)
    assert result["status"] == ("memory_limit" if limit=="memory" else "timeout")
    assert process.killed


@pytest.mark.parametrize("settings", [{"max_pair_evaluations":0}, {"wall_seconds":np.nan}, {"wall_seconds":-1}])
def test_invalid_budgets_reject_before_work(cli, matrix, settings):
    with pytest.raises(ValueError):
        cli.run(matrix.path, matrix.python, matrix.root/"out", **settings)


def test_successful_exit_without_memory_observation_is_not_accepted(cli, matrix, monkeypatch):
    process = SimpleNamespace(pid=123, returncode=0, poll=lambda:0, wait=lambda:0)
    monkeypatch.setattr(cli.subprocess, "Popen", lambda *args, **kwargs:process)
    monkeypatch.setattr(cli.psutil, "Process", lambda pid:SimpleNamespace(children=lambda **kwargs:[]))
    result = cli.run_worker(["never-executed"], cwd=matrix.root, env={}, output=matrix.root/"worker.out", timeout_seconds=120)
    assert result["status"] == "memory_monitor_unavailable"
    assert result["rss_samples"] == 0


def test_historical_windows_data_paths_relocate_without_rewriting_provenance(cli, matrix, monkeypatch):
    clone = matrix.root/"relocated_checkout"
    clone.mkdir()
    historical_root = "Z:\\absent_origin\\pincer-catmech-ai\\data\\"
    payload = deepcopy(matrix.payload)
    original_signatures = [row["signature_sha256"] for row in payload["state_records"]]
    for row in payload["state_records"]:
        for key, relative in (("source_xyz", "structures/Fe_mock_fixture/active/best_found.xyz"),
                              ("source_metadata_path", "structures/Fe_mock_fixture/active/metadata.json")):
            destination = clone/"data"/relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(row[key], destination)
            row[key] = historical_root+relative.replace("/", "\\")
        for item in row["native_evidence"]:
            relative = f"phase3/spin/matrix/multiplicity_{row['multiplicity']}/response.json"
            destination = clone/"data"/relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(item["path"], destination)
            item["path"] = historical_root+relative.replace("/", "\\")
    child = clone/"data/phase3/spin/matrix/summary.json"
    cli.dump(child, payload)
    aggregate = clone/"data/phase3/spin/summary.json"
    cli.dump(aggregate, dict(active_matrix_summary=historical_root+"phase3\\spin\\matrix\\summary.json",
                            active_matrix_summary_sha256=cli.digest(child)))
    before = child.read_bytes()
    monkeypatch.setattr(cli, "REPO", clone)
    ready = cli.inspect_matrix(aggregate, matrix.python, identity_checker=identity_pass)
    assert ready["status"] == "eligible_pair_ready"
    assert ready["selected_seed"]["source_state_signatures"] == original_signatures
    assert child.read_bytes() == before  # Original absolute strings/signature material remain immutable.
    tampered = clone/"data/phase3/spin/matrix/multiplicity_3/response.json"
    tampered.write_text(tampered.read_text()+"changed")
    assert cli.inspect_matrix(aggregate, matrix.python, identity_checker=identity_pass)["status"] == "no_eligible_pair"


def test_existing_literal_evidence_does_not_fall_back_to_a_different_clone_copy(cli, tmp_path, monkeypatch):
    source, clone = tmp_path/"original/data/result.json", tmp_path/"clone"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"original bytes")
    local = clone/"data/result.json"
    local.parent.mkdir(parents=True)
    local.write_bytes(b"different clone bytes")
    monkeypatch.setattr(cli, "REPO", clone)
    assert cli.resolve_evidence_path(source) == source
    assert cli.digest(source) == hashlib.sha256(b"original bytes").hexdigest()


@pytest.mark.parametrize("historical", ["Z:/absent/data/../secret.json", "Z:/absent/data/folder/../../secret.json"])
def test_relocation_rejects_escaping_suffixes(cli, tmp_path, monkeypatch, historical):
    monkeypatch.setattr(cli, "REPO", tmp_path)
    with pytest.raises(ValueError, match="Unsafe historical data path"):
        cli.resolve_evidence_path(historical)


def test_ambiguous_data_suffixes_are_not_guessed(cli, tmp_path, monkeypatch):
    for relative in ("data/parent/data/result.json", "data/result.json"):
        path = tmp_path/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("same bytes")
    monkeypatch.setattr(cli, "REPO", tmp_path)
    with pytest.raises(ValueError, match="Ambiguous historical data path"):
        cli.resolve_evidence_path("Z:/absent/data/parent/data/result.json")


def test_missing_nondata_external_paths_are_not_relocated(cli, tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "REPO", tmp_path)
    local = tmp_path/"data/engine.bin"
    local.parent.mkdir()
    local.write_bytes(b"must not substitute an unrelated runtime")
    with pytest.raises(FileNotFoundError):
        cli.digest("Z:/absent/environment/engine.bin")


def test_runtime_hash_reads_do_not_relocate_even_if_an_environment_is_named_data(cli, tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "REPO", tmp_path)
    local = tmp_path/"data/python.exe"
    local.parent.mkdir()
    local.write_bytes(b"must not substitute the requested runtime")
    with pytest.raises(FileNotFoundError):
        cli.digest("Z:/absent/environment/data/python.exe")


def test_digest_streams_large_inputs_without_read_bytes(cli, tmp_path, monkeypatch):
    payload = b"test chunk boundary\x00"*(200000)
    path = tmp_path/"large.psi-cache"
    path.write_bytes(payload)
    monkeypatch.setattr(Path, "read_bytes", lambda self:pytest.fail("Whole-file reads defeat bounded memory hashing"))
    assert cli.digest(path) == hashlib.sha256(payload).hexdigest()
