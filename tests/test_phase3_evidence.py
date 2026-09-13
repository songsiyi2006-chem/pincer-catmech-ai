"""Reject corrupted, mixed, or ambiguous native evidence before aggregation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import zipfile

import pytest

SPEC = importlib.util.spec_from_file_location("advanced_evidence", Path(__file__).parents[1] / "scripts/run_advanced_campaign.py")
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def fixture_archive(tmp_path, *, multiple):
    archives, members = [], []
    for index in range(2 if multiple else 1):
        name = f"native_evidence_{index+1:03d}.zip" if multiple else "native_evidence.zip"
        payload = f"actual native version {index}".encode()
        with zipfile.ZipFile(tmp_path / name, "w") as archive:
            archive.writestr("native/output", payload)
        archives.append({"name": name, "sha256": RUNNER.sha(tmp_path / name),
                         "bytes": (tmp_path / name).stat().st_size, "member_count": 1})
        members.append({"archive": name, "member": "native/output", "bytes": len(payload),
                        "sha256": hashlib.sha256(payload).hexdigest()})
    document = {"schema": "native_evidence_multiarchive_v1", "archives": archives, "members": members} if multiple else [
        {"path": x["member"], "bytes": x["bytes"], "sha256": x["sha256"]} for x in members]
    (tmp_path / "native_manifest.json").write_text(json.dumps(document), encoding="utf-8")
    return document


@pytest.mark.parametrize("multiple", [False, True])
def test_every_native_member_version_verified(tmp_path, multiple):
    fixture_archive(tmp_path, multiple=multiple)
    audit = RUNNER.verify_archive(tmp_path)
    assert audit["verified"]
    assert audit["members"] == (2 if multiple else 1)


@pytest.mark.parametrize("mutation", ["member_sha", "archive_sha", "duplicate_member", "unknown_archive", "member_count", "extra_member"])
def test_corruption_or_ambiguous_provenance_rejected(tmp_path, mutation):
    document = fixture_archive(tmp_path, multiple=True)
    if mutation == "member_sha":
        document["members"][0]["sha256"] = "0" * 64
    elif mutation == "archive_sha":
        document["archives"][0]["sha256"] = "0" * 64
    elif mutation == "duplicate_member":
        document["members"].append(document["members"][0].copy())
    elif mutation == "unknown_archive":
        document["members"][0]["archive"] = "missing.zip"
    elif mutation == "member_count":
        document["archives"][0]["member_count"] = 19
    else:
        path = tmp_path / document["archives"][0]["name"]
        with zipfile.ZipFile(path, "a") as archive:
            archive.writestr("extra", b"unlisted")
        document["archives"][0].update(sha256=RUNNER.sha(path), bytes=path.stat().st_size)
    (tmp_path / "native_manifest.json").write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(RuntimeError):
        RUNNER.verify_archive(tmp_path)


@pytest.mark.parametrize("fault", [None, "matrix", "not_searched", "invented_minimum", "pair_count"])
def test_target_mecp_gate_remains_separate_from_vertical_scan(fault):
    receipt = dict(completed_utc="fixture", readiness={"matrix_sha256": "bound", "matrix_run_complete": True,
        "eligible_pairs": []}, status="no_eligible_pair", quantum_pair_evaluations=0,
        quantum_state_attempts=0, first_order_crossing=False, minimum_verified=False,
        attempt_result="fixture.json", attempt_result_sha256="fixture")
    if fault == "matrix":
        receipt["readiness"]["matrix_sha256"] = "another"
    elif fault == "not_searched":
        receipt["status"] = "eligible_pair_ready"
    elif fault == "invented_minimum":
        receipt.update(first_order_crossing=True, minimum_verified=True)
    elif fault == "pair_count":
        receipt["quantum_pair_evaluations"] = 1
    if fault is None:
        RUNNER.validate_mecp_receipt(receipt, "bound")
    else:
        with pytest.raises(RuntimeError):
            RUNNER.validate_mecp_receipt(receipt, "bound")


def test_existing_prospective_output_path_is_not_relocated_to_historical_repo(tmp_path):
    path = tmp_path / "data/independent/summary.json"
    path.parent.mkdir(parents=True)
    path.write_text("{}")
    assert RUNNER.repo_artifact(path) == path
