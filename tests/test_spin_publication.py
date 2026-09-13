"""Scientific evidence must remain published when only large caches are omitted."""
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("spin_publication", Path(__file__).parents[1] / "scripts/publish_phase3_spin_evidence.py")
PUBLICATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PUBLICATION)


@pytest.mark.parametrize("name,expected", [("psi.123.97", "PSIF_DFSCF_BJ"),
    ("psi.123.64", "PSIF_LIBDIIS"), ("psi.123.11", None),
    ("psi4.out", None), ("response.json", None), ("psi.result.97", None)])
def test_only_identified_numbered_intermediates_qualify(name, expected):
    assert PUBLICATION.cache_kind(name) == expected


def prepare(tmp_path):
    (tmp_path / "summary.json").write_text(json.dumps({"run_complete": True, "target_states_recorded": 54}))
    (tmp_path / "psi4.out").write_text("Native output fixture\n")
    (tmp_path / "psi.123.97").write_bytes(b"reconstructible cache")
    (tmp_path / "psi.123.11").write_text("gradient fixture\n")
    return PUBLICATION.publish(tmp_path)


def test_cache_local_retention_and_portable_native_verification(tmp_path):
    receipt = prepare(tmp_path)
    assert len(receipt["local_only_caches"]) == 1
    assert len(receipt["published_files"]) == 3
    assert PUBLICATION.verify(tmp_path, local_caches=True)["verified"]
    (tmp_path / "psi.123.97").unlink()  # Simulate a clone lacking omitted cache.
    assert PUBLICATION.verify(tmp_path)["verified"]
    with pytest.raises(FileNotFoundError):
        PUBLICATION.verify(tmp_path, local_caches=True)


def test_missing_scientific_file_or_unlisted_new_file_rejected(tmp_path):
    prepare(tmp_path)
    (tmp_path / "unlisted_gradient").write_text("new evidence")
    with pytest.raises(RuntimeError, match="membership"):
        PUBLICATION.verify(tmp_path)
    (tmp_path / "unlisted_gradient").unlink()
    (tmp_path / "psi4.out").write_text("changed evidence")
    with pytest.raises(RuntimeError, match="SHA/size"):
        PUBLICATION.verify(tmp_path)


def test_incomplete_physical_campaign_cannot_be_frozen(tmp_path):
    (tmp_path / "summary.json").write_text(json.dumps({"run_complete": False, "target_states_recorded": 53}))
    with pytest.raises(RuntimeError, match="Finish"):
        PUBLICATION.publish(tmp_path)


def test_local_verification_rejects_unindexed_new_cache(tmp_path):
    prepare(tmp_path)
    (tmp_path / "psi.456.64").write_bytes(b"new local cache")
    with pytest.raises(RuntimeError, match="cache inventory"):
        PUBLICATION.verify(tmp_path, local_caches=True)
