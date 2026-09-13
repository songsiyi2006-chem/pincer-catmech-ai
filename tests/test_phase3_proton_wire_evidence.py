"""Synthetic byte fixtures test packaging only; no invented QC observations."""
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys
import zipfile

import pytest

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("wire_evidence", REPO / "scripts" / "finalize_phase3_proton_wire.py")
EVIDENCE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EVIDENCE)


def fixture_tree(tmp_path):
    output, scratch = tmp_path / "published", tmp_path / "raw"
    output.mkdir()
    scratch.mkdir()
    rng = random.Random(918)
    for run in ("attempt_a", "attempt_b"):
        (scratch / run).mkdir()
        for index in range(4):
            (scratch / run / f"byte_fixture_{index}.bin").write_bytes(rng.randbytes(240))
    return output, scratch


def test_all_raw_members_preserved_and_archives_split_below_bound(tmp_path):
    output, scratch = fixture_tree(tmp_path)
    manifest = EVIDENCE.archive_native_tree(scratch, output, max_archive_bytes=1000)
    assert len(manifest["archives"]) > 2
    assert all(row["bytes"] < 1000 for row in manifest["archives"])
    actual = {}
    for record in manifest["members"]:
        with zipfile.ZipFile(output / record["archive"]) as archive:
            actual[record["member"]] = archive.read(record["member"])
    expected = {path.relative_to(scratch).as_posix(): path.read_bytes() for path in scratch.rglob("*") if path.is_file()}
    assert actual == expected
    assert EVIDENCE.verify_archives(output)["archive_members_verified"] == 8
    with pytest.raises(FileExistsError, match="immutable"):
        EVIDENCE.archive_native_tree(scratch, output, max_archive_bytes=1000)


def test_corrupted_archive_is_detected_without_touching_raw_evidence(tmp_path):
    output, scratch = fixture_tree(tmp_path)
    manifest = EVIDENCE.archive_native_tree(scratch, output, max_archive_bytes=1000)
    raw_before = {p: EVIDENCE.sha256(p) for p in scratch.rglob("*") if p.is_file()}
    target = output / manifest["archives"][0]["name"]
    payload = bytearray(target.read_bytes())
    payload[45] ^= 1
    target.write_bytes(payload)
    with pytest.raises(ValueError, match="SHA256"):
        EVIDENCE.verify_archives(output)
    assert raw_before == {p: EVIDENCE.sha256(p) for p in raw_before}


def test_finalization_is_idempotent_and_does_not_forge_execution_provenance(tmp_path):
    output, scratch = fixture_tree(tmp_path)
    (output / "summary.json").write_text(json.dumps({"status": "completed_byte_fixture_only", "attempts": [],
        "accepted_TS_count": 0, "activation_free_energy_barrier_count": 0}), encoding="utf-8")
    receipt = EVIDENCE.finalize_evidence(output, scratch, max_archive_bytes=1000)
    assert receipt["new_quantum_calculations_performed"] is False
    assert receipt["executed_driver_source_status"] == "not_recorded_by_this_run"
    assert receipt["completed_path_runs"] == 0 and receipt["actual_band_image_rows"] == 0
    before = {path: EVIDENCE.sha256(path) for path in output.rglob("*") if path.is_file()}
    repeated = EVIDENCE.finalize_evidence(output, scratch, max_archive_bytes=1000)
    assert repeated["already_finalized"] and repeated["existing_receipt_unchanged"]
    assert before == {path: EVIDENCE.sha256(path) for path in before}


def test_output_must_not_be_archived_into_itself(tmp_path):
    with pytest.raises(ValueError, match="non-nested"):
        EVIDENCE.separate_paths(tmp_path / "raw" / "output", tmp_path / "raw")


def test_oversize_single_member_fails_without_a_false_manifest(tmp_path):
    output, scratch = fixture_tree(tmp_path)
    (scratch / "attempt_a" / "too_large.bin").write_bytes(random.Random(8).randbytes(4000))
    with pytest.raises(ValueError, match="Single compressed member"):
        EVIDENCE.archive_native_tree(scratch, output, max_archive_bytes=1000)
    assert (scratch / "attempt_a" / "too_large.bin").stat().st_size == 4000
    assert not (output / "native_manifest.json").exists()


def test_driver_finalization_cli_needs_no_xtb_and_uses_custom_directories(tmp_path):
    output, scratch = fixture_tree(tmp_path)
    (output / "summary.json").write_text(json.dumps({"status": "completed_byte_fixture_only", "attempts": [],
        "accepted_TS_count": 0, "activation_free_energy_barrier_count": 0}), encoding="utf-8")
    result = subprocess.run([sys.executable, str(REPO / "scripts" / "run_phase3_proton_wire.py"),
        "--output", str(output), "--scratch", str(scratch), "--finalize-evidence"],
        cwd=tmp_path, capture_output=True, text=True, encoding="utf-8", timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    receipt = json.loads((output / "verification.json").read_text(encoding="utf-8"))
    assert receipt["archive_members_verified"] == 8
    assert receipt["new_quantum_calculations_performed"] is False
