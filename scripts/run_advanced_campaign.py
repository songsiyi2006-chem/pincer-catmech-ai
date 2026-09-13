"""Run, audit and benchmark Phase 3 without manufacturing chemical labels.

--execute-physical explicitly launches the sequential, bounded native runners.
The default reuses auditable native evidence and reruns the full numerical grid.
Software-fixture TOF/control plots are never presented as catalyst predictions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from pincer_catmech.generators.combinatorial_pincer import BACKBONES, METALS, SUBSTITUENTS
from pincer_catmech.kinetics.master_kinetics import (
    MasterKinetics, MasterThermochemistry, SPECIES, STEP_NAMES, METAL, S, software_fixture_initial,
    software_fixture_rates,
)
from pincer_catmech.kinetics.microkinetics import ComputedFreeEnergy

TEMPERATURES = (360., 370., 380., 383.15, 390., 400., 410., 420.)
BASE_GRID = (.01, .025, .05, .10, .20)
FIXTURE = "software_fixture_not_chemical_prediction"


def sha(path):
    with Path(path).open("rb") as handle:
        return stream_sha(handle)


def stream_sha(handle):
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(1024*1024), b""):
        digest.update(block)
    return digest.hexdigest()


def repo_artifact(path):
    """Resolve historical absolute data paths to this checkout for portability."""
    normalized = str(path).replace("\\", "/")
    candidate = Path(normalized)
    if candidate.is_absolute() and candidate.is_file():
        return candidate
    if "/data/" in normalized:
        return REPO / ("data/"+normalized.split("/data/", 1)[1])
    return candidate if candidate.is_absolute() else REPO / candidate


def dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


def write_csv(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows supplied for {path}")
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def software_benchmark(output):
    """40 batch IVPs and 40 exact tangent IVPs, each with 500 reported times."""
    out = output / "kinetics/software_verification"
    out.mkdir(parents=True, exist_ok=True)
    summaries, trajectories, controls = [], [], []
    arrays = {}
    max_jacobian_error = 0.
    for temperature in TEMPERATURES:
        for base in BASE_GRID:
            model = MasterKinetics(software_fixture_rates(temperature),
                                  tbuoh_activity=.1, allow_software_fixture=True)
            initial = software_fixture_initial(free_base_equiv=base)
            trajectory = model.integrate(initial, end_time=3600., points=500)
            control = model.control_coefficients(trajectory)
            tag = f"T{temperature:g}_base{base:g}"
            arrays[tag+"_concentrations_M"] = trajectory.concentrations
            arrays[tag+"_rate_control"] = control["rate_control"]
            arrays[tag+"_selectivity_control"] = control["instantaneous_selectivity_control"]
            arrays[tag+"_cumulative_product_control"] = control["cumulative_product_control"]
            summaries.append({
                "evidence_kind": FIXTURE, "temperature_K": temperature,
                "initial_free_base_equiv": base, "tbuoh_activity_assumed": .1,
                "points": 500, "time_end_s": 3600.,
                "metal_error_M": trajectory.metal_error_M,
                "base_equivalent_error_M": trajectory.base_equivalent_error_M,
                "metal_sensitivity_residual_M": control["maximum_metal_sensitivity_residual_M"],
                "minimum_concentration_M": trajectory.minimum_concentration_M,
                "terminal_fixture_TOF_s-1": trajectory.tof_per_second[-1],
                "terminal_fixture_product_M": trajectory.concentrations[-1, 14],
                "terminal_fixture_active_metal_fraction": trajectory.active_metal_fraction[-1],
                "steady_state_Campbell_DRC": "not_established",
            })
            # Independent central difference audit covers zero and sampled states.
            for c in (initial, trajectory.concentrations[249], trajectory.concentrations[-1]):
                h = 1e-6
                fd = np.column_stack([(model.rhs(0, c+h*np.eye(18)[j])-model.rhs(0, c-h*np.eye(18)[j]))/(2*h) for j in range(18)])
                max_jacobian_error = max(max_jacobian_error, float(np.max(np.abs(fd-model.jacobian(0, c)))))
            if temperature == 383.15 and base == .05:
                arrays["baseline_concentration_sensitivities_M"] = control["concentration_sensitivities"]
                for q, (t, c, r) in enumerate(zip(trajectory.time, trajectory.concentrations, trajectory.fluxes)):
                    trajectories.append({"evidence_kind": FIXTURE, "time_s": t,
                        **{name+"_M": c[j] for j, name in enumerate(SPECIES)},
                        "fixture_TOF_s-1": trajectory.tof_per_second[q],
                        "active_metal_fraction": trajectory.active_metal_fraction[q],
                        "cleavage_fragment_formed_M": trajectory.cleavage_fragment_formed_M[q]})
                    for j, step in enumerate(STEP_NAMES):
                        row = {"evidence_kind": FIXTURE, "time_s": t, "step": step}
                        for key in ("rate_control", "instantaneous_selectivity_control", "cumulative_product_control"):
                            value = control[key][q, j]
                            row[key] = float(value) if np.isfinite(value) else ""
                        controls.append(row)
            print(json.dumps({"software_grid": tag, "metal_error_M": trajectory.metal_error_M}), flush=True)
    arrays["time_s"] = trajectory.time
    arrays["species"] = np.array(SPECIES)
    arrays["steps"] = np.array(STEP_NAMES)
    arrays["evidence_kind"] = np.array(FIXTURE)
    np.savez_compressed(out / "full_40_grid_500_timepoints.npz", **arrays)
    write_csv(out / "grid_summary.csv", summaries)
    write_csv(out / "baseline_trajectory.csv", trajectories)
    write_csv(out / "baseline_control_coefficients.csv", controls)
    write_csv(out / "stoichiometry.csv", [{"species": name, "metal_count": METAL[i],
              **{step: S[i, j] for j, step in enumerate(STEP_NAMES)}} for i, name in enumerate(SPECIES)])
    result = {"evidence_kind": FIXTURE, "grid_count": len(summaries), "timepoints_per_grid": 500,
        "temperature_K": list(TEMPERATURES), "initial_free_base_equiv": list(BASE_GRID),
        "maximum_jacobian_central_difference_error": max_jacobian_error,
        "maximum_metal_error_M": max(x["metal_error_M"] for x in summaries),
        "maximum_base_equivalent_error_M": max(x["base_equivalent_error_M"] for x in summaries),
        "maximum_metal_sensitivity_residual_M": max(x["metal_sensitivity_residual_M"] for x in summaries),
        "physical_catalyst_predictions": 0,
        "limitations": ["Arbitrary rates test solver behavior only", "Free base is not total added tBuOK",
          "Neutral tBuOH reservoir activity 0.1 is a software fixture assumption",
          "Batch sensitivities are not verified steady-state Campbell DRC",
          "No activation/deactivation barrier inferred from reaction free energies"],
        "files": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "summary.json"}}
    if result["maximum_metal_error_M"] >= 1e-10 or max_jacobian_error >= 1e-5:
        raise RuntimeError("Required numerical accuracy threshold was not met")
    dump(out / "summary.json", result)
    return result


def verify_archive(directory):
    manifest = directory / "native_manifest.json"
    if not manifest.exists():
        raise RuntimeError(f"Missing raw evidence archive/manifest in {directory}")
    document = json.loads(manifest.read_text(encoding="utf-8"))
    if isinstance(document, list):
        archives = [{"name": "native_evidence.zip"}]
        entries = [{"archive": "native_evidence.zip", "member": x["path"],
                    "sha256": x["sha256"], "bytes": x["bytes"]} for x in document]
    elif document.get("schema") == "native_evidence_multiarchive_v1":
        archives, entries = document["archives"], document["members"]
    else:
        raise RuntimeError("Unsupported native evidence manifest schema")
    if not archives or not entries:
        raise RuntimeError("A completed native stage requires a nonempty archive manifest")
    names = [x["name"] for x in archives]
    keys = [(x["archive"], x["member"]) for x in entries]
    if len(names) != len(set(names)) or len(keys) != len(set(keys)):
        raise RuntimeError("Duplicate manifest archive/member identity")
    if set(x[0] for x in keys) != set(names):
        raise RuntimeError("Manifest contains unknown or empty archive")
    audits = []
    for record in archives:
        name = record["name"]
        if Path(name).name != name or "\\" in name:
            raise RuntimeError("Archive name must be a local filename")
        archive = directory / name
        digest = sha(archive)
        if ("sha256" in record and digest != record["sha256"]) or ("bytes" in record and archive.stat().st_size != record["bytes"]):
            raise RuntimeError(f"Native archive SHA/size mismatch: {name}")
        expected = [x for x in entries if x["archive"] == name]
        if "member_count" in record and record["member_count"] != len(expected):
            raise RuntimeError("Native archive member count differs")
        with zipfile.ZipFile(archive) as z:
            if len(z.namelist()) != len(set(z.namelist())):
                raise RuntimeError("Duplicate raw evidence archive members")
            if set(z.namelist()) != {entry["member"] for entry in expected}:
                raise RuntimeError("Raw archive/manifest membership differs")
            # Reading every member to EOF verifies CRC as well as SHA256.
            for entry in expected:
                with z.open(entry["member"]) as handle:
                    member_digest = stream_sha(handle)
                if member_digest != entry["sha256"] or z.getinfo(entry["member"]).file_size != entry["bytes"]:
                    raise RuntimeError(f"Raw evidence SHA/size mismatch: {entry['member']}")
        audits.append({"name": name, "members": len(expected), "sha256": digest, "bytes": archive.stat().st_size})
    return {"members": len(entries), "archives": audits, "manifest_sha256": sha(manifest), "verified": True}


def sibling_script(name):
    spec = importlib.util.spec_from_file_location("phase3_" + name, REPO / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def require_completed_stages(output, summaries):
    """A receipt's existence is insufficient; require the completed stage schema.

    Failed physical searches may be complete. Their negative scientific outcome
    remains separate from successful execution and portable evidence admission.
    """
    mandatory = {
        "solvation": ("finished_utc", "protocol", "clusters", "selected", "association_thermochemistry",
                       "converged_optimizations", "native_archive", "native_manifest"),
        "spin": ("active_matrix_summary", "active_matrix_summary_sha256", "run_complete",
                 "planned_target_states", "target_states_recorded"),
        "proton_wire": ("status", "finished_utc", "method", "solvent", "solvation_state", "attempts",
                        "native_evidence", "accepted_TS_count", "activation_free_energy_barrier_count"),
        "egnn": ("status", "computation_finished", "model_sha256", "dataset_manifest_sha256",
                 "split_sha256", "run_directory", "symmetry"),
        "kinetics_software": ("evidence_kind", "grid_count", "timepoints_per_grid", "temperature_K",
                              "initial_free_base_equiv", "physical_catalyst_predictions", "files"),
    }
    for name, keys in mandatory.items():
        record = summaries.get(name)
        if not isinstance(record, dict) or any(key not in record for key in keys):
            raise RuntimeError(f"Required completed {name} receipt fields are missing")
    solv = summaries["solvation"]
    if not solv["finished_utc"] or (solv["protocol"].get("method"), solv["protocol"].get("solvent"),
                                   solv["protocol"].get("solvation_state")) != ("GFN2-xTB", "toluene", "gsolv"):
        raise RuntimeError("Solvation stage is unfinished or uses a different campaign protocol")
    clusters = solv["clusters"]
    if len(clusters) != 9 or {(x["n_solvent"], x["seed"]) for x in clusters} != {
            (n, seed) for n in (1, 2, 3) for seed in range(3)}:
        raise RuntimeError("Solvation stage requires all nine distinct cluster attempts")
    if any(x.get("protocol") != solv["protocol"] for x in clusters):
        raise RuntimeError("Cluster and solvation-stage protocols differ")
    if solv["converged_optimizations"] != sum(x["native_result"].get("converged") is True for x in clusters):
        raise RuntimeError("Solvation convergence count disagrees with individual records")
    thermal = solv["association_thermochemistry"]
    expected_thermal = {(n, t) for n in (1, 2, 3) for t in (298.15, *TEMPERATURES)}
    if (len(solv["selected"]) != 3 or len(thermal) != len(expected_thermal)
            or {(x["n_solvent"], x["temperature_K"]) for x in thermal} != expected_thermal
            or any(x.get("status") != "verified_cluster_and_reference_minima"
                   or not np.isfinite(x["delta_G_qRRHO_1M_kcal_mol"]) for x in thermal)):
        raise RuntimeError("The complete selected-minimum thermochemistry grid is required")
    for name, local in (("native_archive", "native_evidence.zip"), ("native_manifest", "native_manifest.json")):
        if sha(output / "solvation" / local) != solv[name]["sha256"]:
            raise RuntimeError("Solvation summary is not bound to its native archive/manifest")
    proton = summaries["proton_wire"]
    if (not proton["status"].startswith("completed_") or not proton["finished_utc"] or not proton["attempts"]
            or (proton["method"], proton["solvent"], proton["solvation_state"]) != ("GFN2-xTB", "toluene", "gsolv")):
        raise RuntimeError("Proton-wire stage must be a completed, recorded campaign attempt")
    if any(x.get("status") in (None, "running", "pending") for x in proton["attempts"]):
        raise RuntimeError("Proton-wire attempt remains pending")
    if (proton["accepted_TS_count"] != sum(x.get("accepted_TS") is True for x in proton["attempts"])
            or proton["activation_free_energy_barrier_count"] != sum(x.get("free_energy_barrier_computed") is True for x in proton["attempts"])):
        raise RuntimeError("Proton-wire success counts disagree with individual outcomes")
    spin = summaries["spin"]
    if spin["run_complete"] is not True or spin["planned_target_states"] != 54 or spin["target_states_recorded"] != 54:
        raise RuntimeError("The full 54-slot Fe/Co/Mn spin matrix has not finished")
    egnn = summaries["egnn"]
    if (egnn["status"] != "auxiliary_benchmark_trained" or egnn["computation_finished"] is not True
            or egnn["symmetry"].get("passed") is not True or egnn["symmetry"].get("model_sha256") != egnn["model_sha256"]):
        raise RuntimeError("EGNN completed-training and checkpoint-bound symmetry receipts are required")
    software = summaries["kinetics_software"]
    if (software["evidence_kind"] != FIXTURE or software["grid_count"] != 40 or software["timepoints_per_grid"] != 500
            or software["temperature_K"] != list(TEMPERATURES) or software["initial_free_base_equiv"] != list(BASE_GRID)
            or software["physical_catalyst_predictions"] != 0 or not software["files"]):
        raise RuntimeError("The complete, explicitly nonphysical 40 by 500 software grid is required")
    for key, threshold in (("maximum_jacobian_central_difference_error", 1e-5),
                           ("maximum_metal_error_M", 1e-10), ("maximum_base_equivalent_error_M", 1e-10),
                           ("maximum_metal_sensitivity_residual_M", 1e-10)):
        value = software.get(key)
        if not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0 or value >= threshold:
            raise RuntimeError(f"Software grid numerical acceptance failed: {key}")
    for filename, expected in software["files"].items():
        if sha(output / "kinetics/software_verification" / filename) != expected:
            raise RuntimeError(f"Software grid artifact SHA mismatch: {filename}")


def verify_spin_matrix(spin):
    matrix = spin["matrix"]
    expected = {(f"{m}_{b}_{s}", mult) for m in ("Fe", "Co", "Mn") for b in BACKBONES
                for s in SUBSTITUENTS for mult in (1, 3, 5)}
    records = matrix.get("state_records", [])
    if (matrix.get("run_complete") is not True or matrix.get("planned_target_states") != 54
            or len(records) != 54 or {(x["catalyst_id"], x["multiplicity"]) for x in records} != expected):
        raise RuntimeError("Spin matrix does not contain the complete distinct 54-state design")
    if not matrix.get("method") or not matrix.get("basis") or matrix.get("chemical_state") != "active":
        raise RuntimeError("Spin matrix requires explicit method, basis and active-state identity")
    for row in records:
        if row.get("status") in (None, "running", "pending"):
            raise RuntimeError("A spin-state attempt remains pending")
        if "protocol" in row:
            protocol = row["protocol"]
            actual = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
            if (actual != row.get("protocol_sha256") or protocol.get("method") != matrix["method"]
                    or protocol.get("basis") != matrix["basis"]):
                raise RuntimeError("Spin state protocol/hash differs from its declared matrix")
        elif row.get("status") == "converged":
            raise RuntimeError("A converged spin state lacks protocol provenance")


def physical_kinetics(input_path, output):
    """Archive, bind and solve fully certified conditional models; never fixtures."""
    return sibling_script("phase3_physical_evidence").physical_kinetics(
        input_path, output, REPO, TEMPERATURES, BASE_GRID)


def collect(output, *, require_physical=True):
    summaries = {}
    receipt_hashes = {}
    for name, relative in (("solvation", "solvation/summary.json"), ("spin", "spin/summary.json"),
                           ("proton_wire", "proton_wire/summary.json"),
                           ("egnn", "egnn/summary.json"), ("kinetics_software", "kinetics/software_verification/summary.json")):
        path = output / relative
        if path.is_file():
            summaries[name] = json.loads(path.read_text(encoding="utf-8"))
            receipt_hashes[name] = sha(path)
        elif require_physical:
            raise RuntimeError(f"Required stage receipt missing: {path}")
        else:
            summaries[name] = {"status": "not_yet_completed"}
    if require_physical:
        require_completed_stages(output, summaries)
    if "converged_optimizations" in summaries["solvation"]:
        summaries["solvation_archive_audit"] = verify_archive(output / "solvation")
    proton = summaries["proton_wire"]
    if "native_evidence" in proton:
        summaries["proton_wire_archive_audit"] = verify_archive(output / "proton_wire")
        verification = output / "proton_wire/verification.json"
        proton["verification"] = json.loads(verification.read_text(encoding="utf-8"))
        proton["verification_sha256"] = sha(verification)
        for attempt in proton["attempts"]:
            result = attempt["result"]
            if sha(repo_artifact(result["path"])) != result["sha256"]:
                raise RuntimeError("Proton-wire result hash mismatch")
        if proton["verification"]["status"] != "passed":
            raise RuntimeError("Proton-wire evidence audit has not passed")
        audit = summaries["proton_wire_archive_audit"]
        verified = proton["verification"]
        if require_physical and (verified.get("completed_path_runs") != len(proton["attempts"])
                or verified.get("archive_members_verified") != audit["members"]
                or verified.get("accepted_TS_count") != proton["accepted_TS_count"]
                or verified.get("activation_free_energy_barrier_count") != proton["activation_free_energy_barrier_count"]
                or {(x["name"], x["sha256"]) for x in verified.get("archives_verified", [])}
                    != {(x["name"], x["sha256"]) for x in audit["archives"]}):
            raise RuntimeError("Proton-wire verification is not bound to these attempts and archives")
    spin = summaries["spin"]
    if "active_matrix_summary" in spin:
        path = repo_artifact(spin["active_matrix_summary"])
        if sha(path) != spin["active_matrix_summary_sha256"]:
            raise RuntimeError("Spin matrix summary hash mismatch")
        spin["matrix"] = json.loads(path.read_text(encoding="utf-8"))
        if require_physical:
            verify_spin_matrix(spin)
            summaries["spin_publication_audit"] = sibling_script("publish_phase3_spin_evidence").verify(output / "spin")
    mecp_path = output / "spin/pincer_mecp_from_matrix_001/summary.json"
    if mecp_path.is_file():
        mecp = json.loads(mecp_path.read_text(encoding="utf-8"))
        if require_physical:
            validate_mecp_receipt(mecp, spin["active_matrix_summary_sha256"])
            if sha(repo_artifact(mecp["attempt_result"])) != mecp["attempt_result_sha256"]:
                raise RuntimeError("MECP attempt receipt hash mismatch")
        summaries["mecp"] = mecp
        receipt_hashes["mecp"] = sha(mecp_path)
    elif require_physical:
        raise RuntimeError("Completed target-MECP search or no-eligible-pair receipt is required")
    egnn = summaries["egnn"]
    if "model_sha256" in egnn:
        run = output / "egnn" / egnn["run_directory"].replace("\\", "/")
        for file, expected in ((run / "model.pt", egnn["model_sha256"]),
                               (run / "split.json", egnn["split_sha256"]),
                               (output / "egnn/dataset_manifest.json", egnn["dataset_manifest_sha256"])):
            if sha(file) != expected:
                raise RuntimeError(f"EGNN artifact hash mismatch: {file}")
        if require_physical:
            sibling_script("train_phase3_egnn").verify_benchmark_evidence(output / "egnn", egnn)
    physical_models = sibling_script("phase3_physical_evidence").collect_physical_models(
        output, REPO, TEMPERATURES, BASE_GRID)
    physical_lookup = {(x["catalyst_id"], float(x["temperature_K"]), float(x["total_added_base_equiv"])): x for x in physical_models}
    rows = []
    for metal in METALS:
        for backbone in BACKBONES:
            for substituent in SUBSTITUENTS:
                catalyst = f"{metal}_{backbone}_{substituent}"
                geometry = REPO / f"data/structures/{catalyst}/active/best_found.xyz"
                for temperature in TEMPERATURES:
                    for base in BASE_GRID:
                        physical = physical_lookup.get((catalyst, temperature, base))
                        rows.append({"catalyst_id": catalyst, "metal": metal, "backbone": backbone,
                            "substituent": substituent, "temperature_K": temperature,
                            "total_added_base_equiv": base, "baseline_total_base_equiv": .05,
                            "solvent": "toluene", "active_geometry_available": geometry.is_file(),
                            "spin_requested_multiplicities": "1;3;5" if metal in {"Mn", "Fe", "Co"} else "4d_reference_not_scanned",
                            "kinetic_status": "computed_conditional_model" if physical else "blocked_missing_18_matched_TS_and_offcycle_state_certificates",
                            "free_base_speciation": physical["speciation_assumption"] if physical else "unresolved",
                            "initial_free_base_M": physical["initial_free_base_M"] if physical else "",
                            "tbuoh_activity": physical["tbuoh_activity"] if physical else "unresolved",
                            "physical_TOF_s-1": physical["terminal_TOF_s-1"] if physical else "", "physical_Campbell_DRC": "", "physical_selectivity_control": "see_certified_model_npz" if physical else ""})
    write_csv(output / "catalyst_temperature_base_readiness.csv", rows)
    state = {"schema_version": 3, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "measured_partial_results_with_explicit_missing_chemical_targets",
        "design_count": 24, "temperature_base_grid_per_design": 40,
        "readiness_rows": len(rows), "accepted_physical_18_step_kinetic_predictions": len(physical_models),
        "conditions": {"solvent": "toluene", "implicit_model_for_xtb": "ALPB gsolv", "baseline_temperature_K": 383.15,
           "total_added_tBuOK_equiv": .05, "total_base_scan_equiv": list(BASE_GRID),
           "free_anion_activity_is_total_base": False, "tbuoh_activity_measured": False},
        "prior_campaign_status_sha256": sha(REPO / "data/CAMPAIGN_STATUS.json"),
        "stages": summaries, "stage_receipt_sha256": receipt_hashes,
        "interpretation": "Verified code and real electronic/thermochemical attempts do not establish a complete catalytic mechanism or journal-level validation."}
    dump(output / "PHASE3_STATUS.json", state)
    return state


def validate_mecp_receipt(mecp, matrix_sha):
    """A vertical scan alone does not satisfy the separately executed MECP gate."""
    required = ("completed_utc", "readiness", "status", "quantum_pair_evaluations",
                "quantum_state_attempts", "first_order_crossing", "minimum_verified",
                "attempt_result", "attempt_result_sha256")
    if any(key not in mecp for key in required) or not mecp["completed_utc"]:
        raise RuntimeError("Incomplete target-MECP receipt")
    ready = mecp["readiness"]
    if ready.get("matrix_sha256") != matrix_sha or ready.get("matrix_run_complete") is not True:
        raise RuntimeError("Target-MECP receipt is not bound to the completed spin matrix")
    if mecp["status"] in {"matrix_incomplete", "invalid_matrix_or_orchestration_failure", "eligible_pair_ready"}:
        raise RuntimeError("Target-MECP orchestration is incomplete")
    if mecp["status"] == "no_eligible_pair":
        if (ready.get("eligible_pairs") != [] or mecp["quantum_pair_evaluations"] != 0 or mecp["quantum_state_attempts"] != 0
                or mecp["minimum_verified"] or mecp["first_order_crossing"]):
            raise RuntimeError("No-eligible-pair receipt conflicts with its recorded evidence")
    elif mecp["quantum_state_attempts"] < 1:
        raise RuntimeError("Eligible target-MECP seed was not physically attempted")
    if mecp["minimum_verified"] and not mecp["first_order_crossing"]:
        raise RuntimeError("MECP curvature cannot replace first-order stationarity")


def spin_gap_plot_groups(matrix):
    """Keep finite failed-screen gaps visible and distinct from eligible diagnostics."""
    eligible, rejected = [], []
    for index, record in enumerate(matrix.get("spin_gaps", [])):
        value = record.get("gap_high_minus_low_hartree")
        if record.get("status") != "computed_vertical_gap" or value is None or not np.isfinite(value):
            continue
        passed = (record.get("diagnostic_label_eligible") is True
                  and record.get("spin_contamination_flag") is False
                  and record.get("intended_catalyst_identity_pass") is True)
        (eligible if passed else rejected).append((index, value * 27.211386245988))
    return eligible, rejected


def plots(output, figure_directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none", "svg.hashsalt": "pincer-phase3"})
    figure_directory.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.savefig(figure_directory / (name+".png"), dpi=220, bbox_inches="tight")
        fig.savefig(figure_directory / (name+".pdf"), bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
        svg = figure_directory / (name+".svg")
        fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
        svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")
        plt.close(fig)

    rows = list(csv.DictReader((output / "solvation/association_energies.csv").open(encoding="utf-8")))
    thermal = list(csv.DictReader((output / "solvation/association_thermochemistry.csv").open(encoding="utf-8")))
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    for n in (1, 2, 3):
        values = [float(x["association_energy_kcal_mol"]) for x in rows if int(x["n_solvent"]) == n and x["association_energy_kcal_mol"]]
        axes[0].scatter([n]*len(values), values, color="#087f8c")
        selected = sorted([x for x in thermal if int(x["n_solvent"]) == n], key=lambda x: float(x["temperature_K"]))
        axes[1].plot([float(x["temperature_K"]) for x in selected], [float(x["delta_G_qRRHO_1M_kcal_mol"]) for x in selected], "o-", label=f"{n} tBuOH")
    axes[0].set(xlabel="Explicit tBuOH molecules", ylabel="Association energy (kcal/mol)", xticks=[1,2,3], title="Sampled ALPB model energies")
    axes[1].set(xlabel="Temperature (K)", ylabel="Association free energy (kcal/mol)", title="Certified minima; qRRHO, 1 M")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.axhline(0, color=".6", lw=.8)
    fig.suptitle("Real micro-solvation calculations; no dehydration TS barrier", fontsize=11)
    fig.tight_layout()
    save(fig, "phase3_micro_solvation")

    baseline = list(csv.DictReader((output / "kinetics/software_verification/baseline_trajectory.csv").open(encoding="utf-8")))
    control = list(csv.DictReader((output / "kinetics/software_verification/baseline_control_coefficients.csv").open(encoding="utf-8")))
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.6))
    time = np.array([float(x["time_s"]) for x in baseline])
    axes[0].plot(time, [float(x["fixture_TOF_s-1"]) for x in baseline], color="#087f8c")
    axes[0].set(xlabel="Time (s)", ylabel="Fixture TOF (1/s)", title="Arbitrary rates; 383.15 K")
    axes[1].plot(time, [float(x["active_metal_fraction"]) for x in baseline], color="#17324d")
    axes[1].set(xlabel="Time (s)", ylabel="Active metal fraction", title="Dimer/poisoning code exercise")
    last = sorted([x for x in control if float(x["time_s"]) == 3600.], key=lambda x: abs(float(x["rate_control"])), reverse=True)[:6]
    labels = [f"Step {STEP_NAMES.index(x['step'])+1}" for x in last]
    axes[2].barh(labels[::-1], [float(x["rate_control"]) for x in last][::-1], color="#bd7647")
    axes[2].set(xlabel="Finite-time TS control", title="Software sensitivity ranking")
    fig.suptitle("SOFTWARE FIXTURE ONLY - not catalyst activity or lifetime predictions", color="#a32822", fontsize=11)
    fig.tight_layout()
    save(fig, "phase3_solver_verification")

    egnn_file = output / "egnn/summary.json"
    if egnn_file.exists():
        egnn = json.loads(egnn_file.read_text(encoding="utf-8"))
        run = output / "egnn" / egnn["run_directory"].replace("\\", "/")
        predictions = [json.loads(line) for line in (run / "predictions.jsonl").read_text(encoding="utf-8").splitlines()]
        test = [x for x in predictions if x["split"] == "test"]
        history = json.loads((run / "learning_history.json").read_text(encoding="utf-8"))
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
        target = [x["target_auxiliary_eV"] for x in test]
        predicted = [x["predicted_auxiliary_eV"] for x in test]
        axes[0].scatter(target, predicted, s=15, alpha=.75, color="#087f8c")
        lo, hi = min(target+predicted), max(target+predicted)
        axes[0].plot([lo, hi], [lo, hi], "--", color=".5", lw=1)
        axes[0].set(xlabel="xTB conformer pair energy difference (eV)", ylabel="EGNN pair prediction (eV)", title=f"Held-out family: {len(test)} pairs")
        axes[1].plot([x["epoch"] for x in history], [x["validation"]["mae_eV"] for x in history], "o-", label="EGNN validation MAE")
        axes[1].axhline(egnn["metrics"]["validation"]["equal_reference_energy_baseline_mae_eV"], color="#bd7647", ls="--", label="Zero-difference baseline")
        axes[1].set(xlabel="Epoch (zero-based)", ylabel="Validation MAE (eV)", title="Checkpoint selected by validation")
        axes[1].legend(fontsize=8, frameon=False)
        fig.suptitle("Actual auxiliary learning benchmark; no barrier/spin-gap training", fontsize=11)
        fig.tight_layout()
        save(fig, "phase3_egnn_auxiliary")

    spin_file = output / "spin/summary.json"
    if spin_file.exists():
        spin = json.loads(spin_file.read_text(encoding="utf-8"))
        matrix = json.loads(repo_artifact(spin["active_matrix_summary"]).read_text(encoding="utf-8")) if "active_matrix_summary" in spin else spin.get("matrix", spin)
        records = matrix.get("state_records", [])
        if records:
            fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
            metals = ("Mn", "Fe", "Co")
            converged = [sum(x["catalyst_id"].startswith(m+"_") and x["status"] == "converged" for x in records) for m in metals]
            failed = [sum(x["catalyst_id"].startswith(m+"_") and x["status"] != "converged" for x in records) for m in metals]
            axes[0].bar(metals, converged, label="SCF converged", color="#087f8c")
            axes[0].bar(metals, failed, bottom=converged, label="Failed or missing", color="#bd7647")
            axes[0].set(ylabel="Multiplicity slots", title="Matched protocol vertical diagnostics")
            axes[0].legend(frameon=False)
            eligible, rejected = spin_gap_plot_groups(matrix)
            if eligible or rejected:
                for values, color, marker, label in (
                    (eligible, "#087f8c", "o", "S² / identity screen passed"),
                    (rejected, "#bd7647", "x", "Raw gap: quality screen failed or unavailable"),
                ):
                    if values:
                        axes[1].scatter(*zip(*values), color=color, marker=marker, s=24, label=label)
                axes[1].axhline(0, color=".5", lw=.8)
                axes[1].set(xlabel="Recorded spin pair index", ylabel="Higher-spin minus lower-spin energy (eV)", title="Small-basis gas-phase raw vertical gaps")
                axes[1].legend(frameon=False, fontsize=7)
            else:
                axes[1].text(.5, .5, "No accepted matched spin gaps\nFailures do not establish absent crossings", ha="center", va="center", transform=axes[1].transAxes)
                axes[1].set_axis_off()
            fig.suptitle("DFT capability diagnostics; no validated MECP or solution spin ranking", fontsize=11)
            fig.tight_layout()
            save(fig, "phase3_spin_diagnostics")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", choices=("all", "software", "collect"), default="all")
    p.add_argument("--output", type=Path, default=REPO / "data/phase3")
    p.add_argument("--execute-physical", action="store_true")
    p.add_argument("--psi4-python", type=Path, help="Existing Psi4 interpreter, required for --execute-physical")
    p.add_argument("--xtb", default=os.environ.get("PINCER_XTB") or shutil.which("xtb"), help="Existing native xTB executable for physical stages")
    p.add_argument("--scratch", type=Path, help="Scratch parent for prospective physical execution")
    p.add_argument("--kinetics-input", type=Path, help="Complete source-certified 18-channel thermochemistry JSON")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--allow-incomplete-stages", action="store_true", help="Development only; deployment does not enable this")
    args = p.parse_args()
    args.output = args.output.resolve()
    if args.scratch is not None:
        args.scratch = args.scratch.resolve()
    if args.psi4_python is not None:
        args.psi4_python = args.psi4_python.resolve()
    if args.xtb is not None:
        args.xtb = str(Path(args.xtb).resolve())
    if args.execute_physical and (args.psi4_python is None or not args.psi4_python.is_file()):
        p.error("--execute-physical requires an existing --psi4-python interpreter")
    if args.execute_physical and (args.xtb is None or not Path(args.xtb).is_file()):
        p.error("--execute-physical requires an existing --xtb executable or PINCER_XTB")
    if args.kinetics_input is not None:
        physical_kinetics(args.kinetics_input, args.output)
    if args.execute_physical:
        scratch = args.scratch or REPO / "work/phase3" / args.output.name
        def run_native(script, options):
            subprocess.run([sys.executable, str(REPO / "scripts" / script), *options], cwd=REPO, check=True)
        jobs = [
            ("run_phase3_spin.py", ["--python", str(args.psi4_python),
                "--first-catalyst", "Fe_bipyridine_pnnoh_iPr", "--basis", "sto-3g",
                "--grid-radial", "35", "--grid-spherical", "110", "--scf-algorithm", "soscf",
                "--soscf-start-convergence", "0.0001", "--e-convergence", "1e-8", "--d-convergence", "1e-6",
                "--output-dir", str(args.output / "spin/pincer_vertical_diagnostic_003"),
                "--scratch-dir", str(scratch / "spin"),
                "--aggregate-path", str(args.output / "spin/summary.json")]),
            ("run_phase3_mecp.py", ["--matrix-summary", str(args.output / "spin/pincer_vertical_diagnostic_003/summary.json"),
                "--python", str(args.psi4_python), "--output-dir", str(args.output / "spin/pincer_mecp_from_matrix_001")]),
            ("publish_phase3_spin_evidence.py", ["--directory", str(args.output / "spin")]),
            ("run_solvation_clusters.py", ["--output", str(args.output / "solvation"),
                "--scratch", str(scratch / "solvation"), "--executable", str(args.xtb)]),
        ]
        for script, options in jobs:
            run_native(script, options)
        proton_dir = args.output / "proton_wire"
        proton_options = ["--output", str(proton_dir), "--scratch", str(scratch / "proton_wire")]
        if (proton_dir / "verification.json").is_file():
            run_native("run_phase3_proton_wire.py", [*proton_options, "--verify-evidence"])
        else:
            if not proton_dir.exists():
                run_native("run_phase3_proton_wire.py", [*proton_options, "--xtb", str(args.xtb),
                    "--cluster-source", str(args.output / "solvation/clusters")])
            previous = json.loads((proton_dir / "summary.json").read_text(encoding="utf-8"))
            if (previous.get("attempts") and not previous["attempts"][0]["accepted_TS"]
                    and previous.get("continuation_count", 0) == 0 and previous["deadline_epoch"]-time.time() >= 120):
                run_native("run_phase3_proton_wire.py", [*proton_options, "--xtb", str(args.xtb), "--continue-first", "--neb-steps", "200"])
            run_native("run_phase3_proton_wire.py", [*proton_options, "--finalize-evidence"])
        run_native("train_phase3_egnn.py", ["--output", str(args.output / "egnn")])
    if args.stage in {"all", "software"}:
        software_benchmark(args.output)
    if args.stage in {"all", "collect"}:
        collect(args.output, require_physical=not args.allow_incomplete_stages)
        if not args.no_plots:
            plots(args.output, REPO / "examples/plots")
    print(json.dumps({"phase3_stage": args.stage, "status": "completed_with_explicit_evidence_boundaries"}), flush=True)


if __name__ == "__main__":
    main()
