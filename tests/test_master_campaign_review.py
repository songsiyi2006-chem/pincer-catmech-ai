"""Small manufactured receipts verify admission rules, never chemical evidence."""
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def runner():
    spec = importlib.util.spec_from_file_location("master_campaign_review", Path(__file__).parents[1] / "scripts/run_advanced_campaign.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def complete_schema(runner, monkeypatch):
    """Only tests receipt shape/count policies; native hashing is explicitly stubbed."""
    monkeypatch.setattr(runner, "sha", lambda path: "fixture-digest")
    protocol = {"method": "GFN2-xTB", "solvent": "toluene", "solvation_state": "gsolv"}
    return {
        "solvation": {"finished_utc": "fixture-time", "protocol": protocol,
            "clusters": [{"n_solvent": n, "seed": s, "protocol": protocol, "native_result": {"converged": True}}
                         for n in (1, 2, 3) for s in range(3)], "selected": [{}, {}, {}],
            "association_thermochemistry": [{"n_solvent": n, "temperature_K": t,
                "delta_G_qRRHO_1M_kcal_mol": 1., "status": "verified_cluster_and_reference_minima"}
                for n in (1, 2, 3) for t in (298.15, *runner.TEMPERATURES)],
            "converged_optimizations": 9, "native_archive": {"sha256": "fixture-digest"},
            "native_manifest": {"sha256": "fixture-digest"}},
        "spin": {"active_matrix_summary": "fixture.json", "active_matrix_summary_sha256": "fixture-digest",
                 "run_complete": True, "planned_target_states": 54, "target_states_recorded": 54},
        "proton_wire": {"status": "completed_bounded_exploration", "finished_utc": "fixture-time",
            "method": "GFN2-xTB", "solvent": "toluene", "solvation_state": "gsolv",
            "attempts": [{"status": "failed", "accepted_TS": False, "free_energy_barrier_computed": False}],
            "native_evidence": {"schema": "fixture-only"}, "accepted_TS_count": 0,
            "activation_free_energy_barrier_count": 0},
        "egnn": {"status": "auxiliary_benchmark_trained", "computation_finished": True,
            "model_sha256": "fixture-model", "dataset_manifest_sha256": "fixture-digest", "split_sha256": "fixture-digest",
            "run_directory": "fixture-run", "symmetry": {"passed": True, "model_sha256": "fixture-model"}},
        "kinetics_software": {"evidence_kind": runner.FIXTURE, "grid_count": 40, "timepoints_per_grid": 500,
            "temperature_K": list(runner.TEMPERATURES), "initial_free_base_equiv": list(runner.BASE_GRID),
            "physical_catalyst_predictions": 0, "files": {"grid.npz": "fixture-digest"},
            "maximum_jacobian_central_difference_error": 1e-12, "maximum_metal_error_M": 1e-15,
            "maximum_base_equivalent_error_M": 1e-15, "maximum_metal_sensitivity_residual_M": 1e-15},
    }


def test_completed_negative_search_schema_is_admissible(complete_schema, runner, tmp_path):
    runner.require_completed_stages(tmp_path, complete_schema)
    assert complete_schema["proton_wire"]["accepted_TS_count"] == 0


@pytest.mark.parametrize("stage", ["solvation", "spin", "proton_wire", "egnn", "kinetics_software"])
def test_existing_empty_receipt_cannot_bypass_required_stage(complete_schema, runner, tmp_path, stage):
    complete_schema[stage] = {}
    with pytest.raises(RuntimeError, match="receipt fields"):
        runner.require_completed_stages(tmp_path, complete_schema)


@pytest.mark.parametrize("mutation,reason", [
    ("cluster_duplicate", "nine distinct"), ("thermal_missing", "thermochemistry grid"),
    ("proton_running", "completed, recorded"), ("proton_fake_success", "success counts"),
    ("spin_incomplete", "54-slot"), ("egnn_unbound_symmetry", "checkpoint-bound"),
    ("software_fake_physical", "explicitly nonphysical"), ("software_nan_error", "numerical acceptance"),
])
def test_incomplete_or_contradictory_stage_rejected(complete_schema, runner, tmp_path, mutation, reason):
    if mutation == "cluster_duplicate":
        complete_schema["solvation"]["clusters"][-1] = complete_schema["solvation"]["clusters"][0]
    elif mutation == "thermal_missing":
        complete_schema["solvation"]["association_thermochemistry"].pop()
    elif mutation == "proton_running":
        complete_schema["proton_wire"]["status"] = "running"
    elif mutation == "proton_fake_success":
        complete_schema["proton_wire"]["accepted_TS_count"] = 1
    elif mutation == "spin_incomplete":
        complete_schema["spin"]["target_states_recorded"] = 53
    elif mutation == "egnn_unbound_symmetry":
        complete_schema["egnn"]["symmetry"]["model_sha256"] = "other-model"
    elif mutation == "software_fake_physical":
        complete_schema["kinetics_software"]["physical_catalyst_predictions"] = 1
    else:
        complete_schema["kinetics_software"]["maximum_metal_error_M"] = float("nan")
    with pytest.raises(RuntimeError, match=reason):
        runner.require_completed_stages(tmp_path, complete_schema)


def test_required_collection_rejects_placeholder_files_before_output(runner, tmp_path, monkeypatch):
    for relative in ("solvation", "spin", "proton_wire", "egnn", "kinetics/software_verification"):
        directory = tmp_path / relative
        directory.mkdir(parents=True)
        (directory / "summary.json").write_text("{}", encoding="utf-8")
    monkeypatch.setattr(runner, "dump", lambda *args: pytest.fail("A failed admission wrote status"))
    monkeypatch.setattr(runner, "write_csv", lambda *args: pytest.fail("A failed admission wrote readiness"))
    with pytest.raises(RuntimeError, match="receipt fields"):
        runner.collect(tmp_path)


def test_spin_complete_count_requires_distinct_actual_matrix_slots(runner):
    records = [{"catalyst_id": f"{m}_{b}_{s}", "multiplicity": k, "status": "missing_geometry"}
               for m in ("Fe", "Co", "Mn") for b in runner.BACKBONES for s in runner.SUBSTITUENTS for k in (1, 3, 5)]
    matrix = {"run_complete": True, "planned_target_states": 54, "state_records": records,
              "method": "pbe", "basis": "sto-3g", "chemical_state": "active"}
    runner.verify_spin_matrix({"matrix": matrix})
    matrix["state_records"][-1] = matrix["state_records"][0]
    with pytest.raises(RuntimeError, match="distinct 54-state"):
        runner.verify_spin_matrix({"matrix": matrix})


def test_spin_converged_state_cannot_omit_protocol_hash(runner):
    records = [{"catalyst_id": f"{m}_{b}_{s}", "multiplicity": k, "status": "missing_geometry"}
               for m in ("Fe", "Co", "Mn") for b in runner.BACKBONES for s in runner.SUBSTITUENTS for k in (1, 3, 5)]
    records[0]["status"] = "converged"
    matrix = {"run_complete": True, "planned_target_states": 54, "state_records": records,
              "method": "pbe", "basis": "sto-3g", "chemical_state": "active"}
    with pytest.raises(RuntimeError, match="lacks protocol"):
        runner.verify_spin_matrix({"matrix": matrix})


def test_empty_multiarchive_is_not_completed_native_evidence(runner, tmp_path):
    (tmp_path / "native_manifest.json").write_text(json.dumps({"schema": "native_evidence_multiarchive_v1",
                                                             "archives": [], "members": []}))
    with pytest.raises(RuntimeError, match="nonempty archive"):
        runner.verify_archive(tmp_path)


def test_spin_plot_retains_but_qualifies_failed_screen_gaps(runner):
    accepted = {"status": "computed_vertical_gap", "gap_high_minus_low_hartree": .1,
                "diagnostic_label_eligible": True, "spin_contamination_flag": False,
                "intended_catalyst_identity_pass": True}
    rejected = {**accepted, "spin_contamination_flag": True, "diagnostic_label_eligible": False}
    absent = {"status": "missing_or_failed_state", "gap_high_minus_low_hartree": None}
    eligible, failed = runner.spin_gap_plot_groups({"spin_gaps": [accepted, rejected, absent]})
    assert eligible == [(0, pytest.approx(2.7211386245988))]
    assert failed == [(1, pytest.approx(2.7211386245988))]
    assert all(value != 0 for _, value in eligible + failed)


def test_readiness_rows_use_validated_conditional_speciation_values(runner, tmp_path, monkeypatch):
    physical = {"catalyst_id": "Fe_bipyridine_pnnoh_iPr", "temperature_K": 383.15,
                "total_added_base_equiv": .05, "initial_free_base_M": .002,
                "tbuoh_activity": .37, "speciation_assumption": "explicit test-only assumption",
                "terminal_TOF_s-1": .004}
    monkeypatch.setattr(runner, "sibling_script", lambda name: SimpleNamespace(collect_physical_models=lambda *args: [physical]))
    monkeypatch.setattr(runner, "sha", lambda path: "fixture-digest")
    saved = []
    monkeypatch.setattr(runner, "write_csv", lambda path, rows: saved.extend(rows))
    monkeypatch.setattr(runner, "dump", lambda *args: None)
    runner.collect(tmp_path, require_physical=False)
    row = next(r for r in saved if r["catalyst_id"] == physical["catalyst_id"]
               and r["temperature_K"] == 383.15 and r["total_added_base_equiv"] == .05)
    assert row["initial_free_base_M"] == .002 and row["tbuoh_activity"] == .37
    assert row["free_base_speciation"] == physical["speciation_assumption"]


def test_physical_wrapper_passes_campaign_catalog_and_grid(runner, tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(runner, "sibling_script", lambda name: SimpleNamespace(
        physical_kinetics=lambda *args: seen.append(args) or []))
    assert runner.physical_kinetics(tmp_path / "input.json", tmp_path) == []
    assert seen == [(tmp_path / "input.json", tmp_path, runner.REPO, runner.TEMPERATURES, runner.BASE_GRID)]
