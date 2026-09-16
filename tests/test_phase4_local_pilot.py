"""Software-only boundary tests. No test datum is a chemical calculation."""
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def cli():
    spec = importlib.util.spec_from_file_location("phase4_local", Path(__file__).resolve().parents[1] / "scripts/run_phase4_local_pilot.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture_response(cli):
    request = dict(**cli.PROTOCOL, symbols=["H", "H"], positions_bohr=[[0., 0., 0.], [0., 0., 1.4]],
                   charge=0, multiplicities=[1], basis="sto-3g")
    response = dict(symbols=request["symbols"], positions_bohr=request["positions_bohr"], charge=0,
                    multiplicities=[1], basis="sto-3g", method="pbe", coordinate_units="bohr",
                    energy_units="hartree", gradient_units="hartree/bohr",
                    geometry_sha256=cli.geometry_digest(request["symbols"], request["positions_bohr"]),
                    scf_e_convergence=1e-8, scf_d_convergence=1e-6,
                    **{key: request[key] for key in ("grid_radial", "grid_spherical", "scf_algorithm", "scf_max_iterations", "soscf_start_convergence")},
                    states=[dict(multiplicity=1, energy_hartree=-1., spin_squared=0., expected_spin_squared=0.,
                                 gradient_hartree_bohr=[[0., 0., 0.], [0., 0., 0.]], scf_converged=True,
                                 reference="rks", gradient_kind="analytic")])
    return request, response


def test_existing_evidence_rejected(cli, tmp_path):
    existing = tmp_path / "prior"
    existing.mkdir()
    (existing / "evidence").write_text("immutable")
    with pytest.raises(FileExistsError):
        cli.paths_for_run(existing, tmp_path / "scratch")
    assert (existing / "evidence").read_text() == "immutable"


@pytest.mark.parametrize("reverse", [False, True])
def test_nested_paths_rejected(cli, tmp_path, reverse):
    a, b = tmp_path / "parent", tmp_path / "parent/child"
    with pytest.raises(ValueError, match="separate"):
        cli.paths_for_run(*(b, a) if reverse else (a, b))


def test_phase3_path_rejected(cli):
    with pytest.raises(ValueError, match="Phase3"):
        cli.paths_for_run(cli.REPO / "data/phase3/new_should_never_exist", cli.REPO / "work/new_scratch")


def test_valid_software_response(cli):
    request, response = fixture_response(cli)
    assert cli.validate_response(response, request)["spin_quality_pass"]


@pytest.mark.parametrize("key,value", [("basis", "def2-svp"), ("geometry_sha256", "wrong"),
    ("charge", 1), ("energy_units", "eV"), ("scf_e_convergence", 1e-6), ("grid_radial", 20)])
def test_wrong_response_identity_protocol_rejected(cli, key, value):
    request, response = fixture_response(cli)
    response[key] = value
    with pytest.raises(ValueError):
        cli.validate_response(response, request)


@pytest.mark.parametrize("key,value", [("scf_converged", False), ("energy_hartree", float("nan")),
    ("spin_squared", float("inf")), ("expected_spin_squared", 2.), ("reference", "uks"),
    ("gradient_hartree_bohr", [[0., 0., 0.]])])
def test_invalid_worker_state_rejected(cli, key, value):
    request, response = fixture_response(cli)
    response["states"][0][key] = value
    with pytest.raises(ValueError):
        cli.validate_response(response, request)


def test_contamination_does_not_become_accepted_label(cli):
    request, response = fixture_response(cli)
    response["states"][0]["spin_squared"] = .2
    assert not cli.validate_response(response, request)["spin_quality_pass"]


def test_native_archive_retains_binary_and_copies_science(cli, tmp_path):
    scratch, output = tmp_path / "scratch", tmp_path / "output"
    scratch.mkdir()
    output.mkdir()
    (scratch / "psi.123.97").write_bytes(b"local temporary integral fixture")
    (scratch / "psi4.out").write_text("mock scientific output; not computed chemistry")
    manifest = cli.archive_native(scratch, output)
    assert manifest["files"] == 2 and manifest["scientific_files_copied"] == 1
    assert (scratch / "psi.123.97").exists()
    assert (output / "native/psi4.out").read_bytes() == (scratch / "psi4.out").read_bytes()
    records = json.loads((output / "native_manifest.json").read_text())["files"]
    assert all(row["sha256"] == cli.digest(scratch / row["relative_path"]) for row in records)


def test_time_limit_never_loosens_scientific_threshold(cli):
    assert sum(row[3] for row in cli.JOBS) == 1140
    assert cli.PROTOCOL["e_convergence"] == 1e-8
    assert cli.PROTOCOL["d_convergence"] == 1e-6
    assert cli.PROTOCOL["threads"] == 2
    assert cli.PROTOCOL["memory_mib"] == 500


def test_no_quantum_launch_after_budget(cli, tmp_path, monkeypatch):
    output, scratch = tmp_path / "out", tmp_path / "scratch"
    output.mkdir()
    scratch.mkdir()
    monkeypatch.setattr(cli.subprocess, "Popen", lambda *a, **k: pytest.fail("No process should launch"))
    args = SimpleNamespace(python=tmp_path / "python.exe", rss_limit_mib=1536)
    result = cli.run_job(("budget", "sto-3g", 1, 240), args, ["H", "H"],
                         np.array([[0., 0., 0.], [0., 0., 1.4]]), output, scratch, cli.time.monotonic() - 1)
    assert result["status"] == "not_started_budget"
    assert result["accepted_chemical_label"] is False
    assert (output / "budget/result.json").exists()


def test_host_memory_gate_preserves_not_started_record(cli, tmp_path, monkeypatch):
    output, scratch = tmp_path / "out", tmp_path / "scratch"
    output.mkdir()
    scratch.mkdir()
    monkeypatch.setattr(cli.subprocess, "Popen", lambda *a, **k: pytest.fail("No process should launch"))
    monkeypatch.setattr(cli.psutil, "virtual_memory", lambda: SimpleNamespace(available=1024**3))
    args = SimpleNamespace(python=tmp_path / "python.exe", rss_limit_mib=1536)
    result = cli.run_job(("memory", "sto-3g", 1, 240), args, ["H", "H"],
                         np.array([[0., 0., 0.], [0., 0., 1.4]]), output, scratch, cli.time.monotonic() + 300)
    assert result["status"] == "not_started_memory"
    assert result["accepted_chemical_label"] is False


@pytest.mark.parametrize("bad_response", [False, True])
def test_mock_worker_evidence_and_nonchemical_acceptance(cli, tmp_path, monkeypatch, bad_response):
    output, scratch = tmp_path / "out", tmp_path / "scratch"
    output.mkdir()
    scratch.mkdir()
    request, response = fixture_response(cli)
    if bad_response:
        response["basis"] = "wrong"
    class MockProcess:
        pid = 123456
        def poll(self):
            return 0
        def wait(self):
            return 0
    def fake_launch(command, **kwargs):
        Path(command[command.index("--result") + 1]).write_text(json.dumps(response))
        kwargs["stdout"].write("MOCK SOFTWARE FIXTURE; NO QUANTUM CALCULATION\n")
        return MockProcess()
    monkeypatch.setattr(cli.subprocess, "Popen", fake_launch)
    monkeypatch.setattr(cli.psutil, "Process", lambda pid: object())
    monkeypatch.setattr(cli.psutil, "virtual_memory", lambda: SimpleNamespace(available=4 * 1024**3))
    args = SimpleNamespace(python=tmp_path / "python.exe", rss_limit_mib=1536)
    result = cli.run_job(("mock", "sto-3g", 1, 240), args, request["symbols"],
                         np.array(request["positions_bohr"]), output, scratch, cli.time.monotonic() + 300)
    assert result["status"] == ("driver_failed" if bad_response else "converged")
    assert result["accepted_chemical_label"] is False
    assert result["response_sha256"] == cli.digest(output / "mock/response.json")
    assert result["worker_stdout_sha256"] == cli.digest(output / "mock/worker.out")
