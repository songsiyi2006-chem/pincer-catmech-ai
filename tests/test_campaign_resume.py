"""Resume-provenance regressions; no quantum program or worker is launched."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

pytest.importorskip("psutil")
SPEC = importlib.util.spec_from_file_location(
    "pincer_conformer_resume_test", Path(__file__).resolve().parents[1]/"scripts/run_conformer_campaign.py")
DRIVER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DRIVER
SPEC.loader.exec_module(DRIVER)


def control(solvent=None):
    return {"status": "finished", "pid": 0, "started_epoch": 1000.0,
            "started_at_utc": "1970-01-01T00:16:40+00:00", "deadline_epoch": 1100.0,
            "requested_hours": 100/3600, "workers": 2, "threads_per_worker": 2,
            **DRIVER.campaign_protocol(solvent)}


@pytest.mark.parametrize("old_solvent,new_solvent", [(None, "toluene"), ("toluene", None)])
def test_resume_cannot_mix_gas_and_toluene_energy_records(tmp_path, monkeypatch, old_solvent, new_solvent):
    monkeypatch.setattr(DRIVER, "REPO", tmp_path)
    folder = tmp_path/"data"/"prior_run"
    folder.mkdir(parents=True)
    old = control(old_solvent)
    original = json.dumps(old).encode()
    (folder/"control.json").write_bytes(original)
    # An explicit rejection fixture. It must never be read as a quantum record.
    evidence = b'{"analytical_rejection_fixture": true, "energy_eV": -1.0}\n'
    (folder/"records.jsonl").write_bytes(evidence)
    monkeypatch.setattr(DRIVER, "ProcessPoolExecutor", lambda **kwargs: pytest.fail("worker creation before protocol gate"))
    with pytest.raises(ValueError, match="scientific protocol"):
        DRIVER.run_campaign(3, 1, 1, "must-not-execute-xtb", resume=True,
                            solvent=new_solvent, campaign_name="prior_run")
    assert (folder/"control.json").read_bytes() == original
    assert (folder/"records.jsonl").read_bytes() == evidence
    assert not (folder/"records_with_identity.jsonl").exists()
    assert not (folder/"datasets").exists()


@pytest.mark.parametrize("field,value", [("method", "GFN1-xTB"), ("solvation_state", "bar1mol"),
                                         ("electronic_temperature_K", 1000.0), ("method", None)])
def test_resume_requires_all_fixed_protocol_fields(field, value):
    previous = control("toluene")
    previous[field] = value
    with pytest.raises(ValueError, match=field):
        DRIVER.validate_resume_protocol(previous, "toluene")


@pytest.mark.parametrize("contents", [None, "{}"])
def test_resume_requires_an_existing_complete_control(tmp_path, monkeypatch, contents):
    monkeypatch.setattr(DRIVER, "REPO", tmp_path)
    folder = tmp_path/"data"/"prior_run"
    folder.mkdir(parents=True)
    if contents is not None:
        (folder/"control.json").write_text(contents)
    monkeypatch.setattr(DRIVER, "ProcessPoolExecutor", lambda **kwargs: pytest.fail("worker creation without protocol"))
    with pytest.raises(ValueError, match="protocol"):
        DRIVER.run_campaign(3, 1, 1, "must-not-execute-xtb", resume=True,
                            solvent="toluene", campaign_name="prior_run")
    assert not (folder/"records_with_identity.jsonl").exists()


def test_same_protocol_resource_changes_preserve_original_deadline(tmp_path, monkeypatch):
    monkeypatch.setattr(DRIVER, "REPO", tmp_path)
    monkeypatch.setattr(DRIVER, "enumerate_specs", lambda: [])
    folder = tmp_path/"data"/"prior_run"
    folder.mkdir(parents=True)
    old = control("toluene")
    (folder/"control.json").write_text(json.dumps(old))

    class NoQuantumPool:
        def __init__(self, **kwargs):
            assert kwargs == {"max_workers": 1}

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, *args, **kwargs):
            pytest.fail("An expired resumed deadline must not submit quantum work")

    monkeypatch.setattr(DRIVER, "ProcessPoolExecutor", NoQuantumPool)
    DRIVER.run_campaign(9, 1, 4, "must-not-execute-xtb", resume=True, solvent="toluene",
                        deadline_epoch=999999999999.0, campaign_name="prior_run")
    result = json.loads((folder/"control.json").read_text())
    assert result["deadline_epoch"] == old["deadline_epoch"]
    assert result["started_epoch"] == old["started_epoch"]
    assert result["requested_hours"] == old["requested_hours"]
    assert result["workers"] == 1
    assert result["threads_per_worker"] == 4
    assert {key: result[key] for key in DRIVER.campaign_protocol("toluene")} == DRIVER.campaign_protocol("toluene")
