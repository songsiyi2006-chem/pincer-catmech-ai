"""Manufactured artifact bundles exercise provenance checks, never chemistry."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

import pytest
import torch


@pytest.fixture
def trainer():
    source = Path(__file__).resolve().parents[1] / "scripts/train_phase3_egnn.py"
    spec = importlib.util.spec_from_file_location("egnn_evidence_review", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def bundle(tmp_path, trainer):
    """A tiny fake trained-artifact container; no optimizer or physical labels."""
    for name in ("structures.jsonl", "pairs.jsonl", "rejected_sources.json", "reference_selection.json"):
        (tmp_path / name).write_text("[]\n", encoding="utf-8")
    manifest = {
        "eligible_structures": 0, "rejected_records": 0, "auxiliary_pairs": 0,
        "composition_groups": 0, "families": {},
        "dataset_files_sha256": {p.name: trainer.digest(p.read_bytes()) for p in tmp_path.iterdir()},
    }
    trainer.atomic_json(tmp_path / "dataset_manifest.json", manifest)
    run = tmp_path / "runs/fixture"
    run.mkdir(parents=True)
    trainer.atomic_json(run / "split.json", {"train": {"indices": [], "family_ids": []}})
    configuration = {"node_feature_dim": 3, "hidden_dim": 8, "num_layers": 1,
                     "heads": list(trainer.DEFAULT_HEADS)}
    model = trainer.PincerEGNN(**configuration)
    model.mark_heads_trained([trainer.AUXILIARY], pair_only=True)
    normalizer = {"mean": [0., 0., .12], "scale": [1., 1., .7],
                  "supported": [False, False, True], "counts": [0, 0, 4]}
    checkpoint = {"configuration": configuration, "model_state_dict": model.state_dict(),
                  "dataset_manifest_sha256": trainer.digest((tmp_path / "dataset_manifest.json").read_bytes())}
    checkpoint.update({"normalization_" + key: torch.tensor(value) for key, value in normalizer.items()})
    torch.save(checkpoint, run / "model.pt")
    report = {"status": "auxiliary_benchmark_trained", "run_directory": "runs/fixture",
              "configuration": configuration, "normalization": normalizer,
              "heads": {head: {"trained": i == 2, "pair_only": i == 2}
                        for i, head in enumerate(trainer.DEFAULT_HEADS)},
              "dataset_manifest_sha256": checkpoint["dataset_manifest_sha256"],
              "model_sha256": trainer.digest((run / "model.pt").read_bytes()),
              "split_sha256": trainer.digest((run / "split.json").read_bytes())}
    trainer.atomic_json(tmp_path / "readiness.json", report)
    return tmp_path, report, checkpoint


def test_valid_evidence_bundle_is_read_only_and_allows_current_source_changes(bundle, trainer, monkeypatch):
    root, report, _ = bundle
    before = {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    monkeypatch.setattr(trainer, "source_fingerprints", lambda: {"current.py": "new-verifier-hash"})
    manifest, checkpoint, split = trainer.verify_benchmark_evidence(root, report)
    assert checkpoint["dataset_manifest_sha256"] == report["dataset_manifest_sha256"]
    assert "structures.jsonl" in manifest["dataset_files_sha256"] and "train" in split
    assert before == {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()}


@pytest.mark.parametrize("relative,reason", [
    ("dataset_manifest.json", "manifest SHA"),
    ("structures.jsonl", "dataset file SHA"),
    ("pairs.jsonl", "dataset file SHA"),
    ("runs/fixture/split.json", "split SHA"),
    ("runs/fixture/model.pt", "checkpoint SHA"),
])
def test_refresh_rejects_changed_artifacts_before_writing_any_summary(bundle, trainer, relative, reason):
    root, report, _ = bundle
    prior_summary = {**report, "review_note": "must remain unchanged after failure"}
    trainer.atomic_json(root / "summary.json", prior_summary)
    before = (root / "summary.json").read_bytes()
    path = root / relative
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match=reason):
        trainer.summarize_benchmark(root)
    assert (root / "summary.json").read_bytes() == before
    assert not (root / "symmetry_verification.json").exists()


@pytest.mark.parametrize("field,reason", [
    ("dataset_manifest_sha256", "different training dataset"),
    ("configuration", "configuration"),
    ("normalization_mean", "normalization mean"),
    ("trained_heads", "trained_heads"),
    ("pair_only_heads", "pair_only_heads"),
])
def test_checkpoint_internal_binding_cannot_be_bypassed_by_updating_outer_hash(bundle, trainer, field, reason):
    root, report, checkpoint = bundle
    changed = copy.deepcopy(checkpoint)
    if field == "dataset_manifest_sha256":
        changed[field] = "0" * 64
    elif field == "configuration":
        changed[field]["num_layers"] = 2
    elif field == "normalization_mean":
        changed[field][2] += 1
    else:
        changed["model_state_dict"][field][0] = True
    path = root / "runs/fixture/model.pt"
    torch.save(changed, path)
    report["model_sha256"] = trainer.digest(path.read_bytes())
    with pytest.raises(ValueError, match=reason):
        trainer.verify_benchmark_evidence(root, report)


def test_legacy_training_fingerprints_survive_verification_with_new_code(bundle, trainer):
    root, report, checkpoint = bundle
    historical = {"scripts/train_phase3_egnn.py": "old-training-code-digest"}
    trainer.atomic_json(root / "summary.json", {**report, "source_files_sha256": historical,
                                               "finished_utc": "original-finish-time"})
    sources, origin, previous = trainer.training_source_provenance(root, report, checkpoint)
    assert sources == historical
    assert "retained" in origin
    assert previous["finished_utc"] == "original-finish-time"


def test_unrelated_prior_summary_cannot_supply_training_source_provenance(bundle, trainer):
    root, report, checkpoint = bundle
    trainer.atomic_json(root / "summary.json", {**report, "model_sha256": "another-model",
                                               "source_files_sha256": {"source.py": "wrong-run"}})
    sources, origin, previous = trainer.training_source_provenance(root, report, checkpoint)
    assert sources is None and previous == {}
    assert "unavailable" in origin


def test_new_checkpoint_training_source_attestation_must_match_report(bundle, trainer):
    root, report, checkpoint = bundle
    checkpoint["training_source_files_sha256"] = {"source.py": "training-source"}
    path = root / "runs/fixture/model.pt"
    torch.save(checkpoint, path)
    report["model_sha256"] = trainer.digest(path.read_bytes())
    with pytest.raises(ValueError, match="training-source provenance"):
        trainer.verify_benchmark_evidence(root, report)
    report["training_source_files_sha256"] = checkpoint["training_source_files_sha256"]
    trainer.verify_benchmark_evidence(root, report)
    sources, origin, _ = trainer.training_source_provenance(root, report, checkpoint)
    assert sources == checkpoint["training_source_files_sha256"]
    assert "embedded in checkpoint" in origin


def test_verify_only_cli_never_prepares_data_or_trains(trainer, monkeypatch, tmp_path, capsys):
    output = tmp_path / "retained benchmark"
    calls = []

    def forbidden(*args, **kwargs):
        pytest.fail("Verification-only CLI must never prepare data or train")

    def summarize(path):
        calls.append(path)
        return {"status": "verified_test_fixture"}

    monkeypatch.setattr(trainer, "build_dataset", forbidden)
    monkeypatch.setattr(trainer, "run_training", forbidden)
    monkeypatch.setattr(trainer, "summarize_benchmark", summarize)
    monkeypatch.setattr(sys, "argv", ["trainer", "--verify-only", "--output", str(output),
                                    "--records", str(tmp_path / "missing-source.jsonl")])
    trainer.main()
    assert calls == [output]
    assert json.loads(capsys.readouterr().out) == {"status": "verified_test_fixture"}
    assert not output.exists()


def test_verify_only_cli_propagates_evidence_failure_without_training(bundle, trainer, monkeypatch):
    output, _, _ = bundle
    changed = output / "pairs.jsonl"
    changed.write_bytes(changed.read_bytes() + b" ")
    before = {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}

    def forbidden(*args, **kwargs):
        pytest.fail("Evidence rejection must not trigger data preparation or training")

    monkeypatch.setattr(trainer, "build_dataset", forbidden)
    monkeypatch.setattr(trainer, "run_training", forbidden)
    monkeypatch.setattr(sys, "argv", ["trainer", "--verify-only", "--output", str(output)])
    with pytest.raises(ValueError, match="dataset file SHA"):
        trainer.main()
    assert before == {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}


def test_prepare_only_and_verify_only_are_mutually_exclusive(trainer, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["trainer", "--prepare-only", "--verify-only"])
    with pytest.raises(SystemExit) as error:
        trainer.main()
    assert error.value.code == 2
    assert "not allowed with argument --prepare-only" in capsys.readouterr().err
