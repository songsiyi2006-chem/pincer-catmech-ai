"""Freeze actual bounded-campaign data after its original compute deadline.

This does not turn failed paths into barriers or incomplete mechanisms into
kinetics. Publication status is separate from scientific outcome acceptance.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from ase.io import read
import psutil
from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now
from pincer_catmech.kinetics.microkinetics import campaign_readiness
from pincer_catmech.models.surrogate_regressor import training_readiness


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record_count(path):
    if not path.is_file():
        return 0
    return sum(bool(line.strip()) for line in path.open(encoding="utf-8"))


def verify_thermochemistry_snapshot():
    """Freeze the exact completed certificates that passed the full exporter."""
    folder = REPO / "data/datasets/thermochemistry"
    protocol = load(folder / "export_protocol.json")
    for name, expected in protocol["files"].items():
        if digest(folder / name) != expected:
            raise ValueError(f"Thermochemistry export changed since validation: {name}")
    selected = load(folder / "selected_certificates.json")
    for key, entry in selected.items():
        certificate = REPO / entry["certificate"]
        if not entry.get("eligible") or not entry.get("chemical_identity", {}).get("retained"):
            raise ValueError(f"Unverified identity in selected certificate: {key}")
        for file_key, hash_key in (("certificate", "certificate_sha256"),
                                   ("accepted_geometry", "accepted_geometry_sha256"),
                                   ("spectrum", "spectrum_sha256")):
            if digest(REPO / entry[file_key]) != entry[hash_key]:
                raise ValueError(f"Selected thermochemistry evidence changed: {key}/{file_key}")
        source = load(certificate)
        if not source.get("finished_utc") or source.get("status") != "verified_minimum" or source.get("accepted") is not True:
            raise ValueError(f"Selected thermochemistry certificate is incomplete: {key}")
    counts = {"selected_catalyst_states": sum(entry["kind"] == "catalyst" for entry in selected.values()),
              "selected_active_catalysts": sum(entry["kind"] == "catalyst" and key.endswith("|active") for key, entry in selected.items()),
              "selected_reference_species": sum(entry["kind"] == "reference" for entry in selected.values())}
    if any(protocol[name] != count for name, count in counts.items()):
        raise ValueError("Thermochemistry selection counts disagree with completed certificates")
    return protocol


def freeze_structures():
    best_path = REPO / "data/datasets/best_conformers.json"
    best = load(best_path)
    portable = {}
    for key, record in best.items():
        metadata = record["metadata"]
        if (key != f"{metadata['catalyst_id']}|{metadata['state']}" or
                metadata["catalyst_id"] != record["catalyst_id"] or metadata["state"] != record["state"]):
            raise ValueError(f"Record identity metadata disagrees: {key}")
        if not record["native_result"].get("converged") or record.get("solvent") != "toluene" or record.get("solvation_state") != "gsolv":
            raise ValueError(f"Missing converged ALPB toluene gsolv provenance: {key}")
        folder = REPO / "data/structures" / record["catalyst_id"] / record["state"]
        target = folder / "best_found.xyz"
        candidates = [target]
        if record.get("geometry"):
            candidates.append(REPO / record["geometry"])
        candidates.append(Path(record["native_result"]["output_directory"]) / "xtbopt.xyz")
        source = next((candidate for candidate in candidates if candidate.is_file()
                       and digest(candidate) == record["optimized_geometry_sha256"]), None)
        if source is None:
            raise ValueError(f"No source geometry with matching SHA256: {key}")
        atoms = read(source)
        if not catalyst_identity_screen(atoms, metadata)["retained"]:
            raise ValueError(f"Intended identity lost: {key}")
        distances = atoms.get_distances(metadata["metal_index"], metadata["donor_indices"]).tolist()
        if not all(1.4 < d < 3.1 for d in distances):
            raise ValueError(f"Donor screen failed: {key}")
        folder.mkdir(parents=True, exist_ok=True)
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        if digest(target) != digest(source):
            raise ValueError("Copied geometry changed")
        atomic_json(folder / "metadata.json", metadata)
        portable[key] = {**record, "geometry": target.relative_to(REPO).as_posix(),
                         "geometry_sha256": digest(target), "metal_donor_distances_A": distances,
                         "metadata_file": (folder / "metadata.json").relative_to(REPO).as_posix(),
                         "metadata_sha256": digest(folder / "metadata.json"),
                         "source_record": "data/datasets/best_conformers.json",
                         "source_record_key": key,
                         "geometry_selection": "lowest electronic energy among sampled identity-preserving optimizations; no Hessian certification implied"}
    atomic_json(REPO / "data/datasets/best_conformers_portable.json", portable)
    descriptors = REPO / "data/datasets/catalyst_descriptors.csv"
    with descriptors.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = list(reader.fieldnames), list(reader)
    fields += [field for field in ("geometry_sha256", "selection_source") if field not in fields]
    for row in rows:
        key = f"{row['catalyst_id']}|active"
        if key in portable:
            row["geometry_path"] = portable[key]["geometry"]
            row["geometry_sha256"] = portable[key]["optimized_geometry_sha256"]
            row["selection_source"] = "data/datasets/best_conformers_portable.json"
        else:
            row.update(geometry_path="", geometry_sha256="", selection_source="")
    with descriptors.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return best, rows


def main():
    original = load(REPO / "data/campaign/control.json")
    alpb = load(REPO / "data/campaign_alpb_toluene/control.json")
    if time.time() < original["deadline_epoch"] or alpb["status"] != "finished":
        raise RuntimeError("Cannot freeze before original deadline and main worker completion")
    tracked = {"run_conformer_campaign.py", "run_neb_campaign.py", "run_minimum_thermochemistry.py",
               "run_reference_thermochemistry.py", "run_condensation_search.py", "monitor_campaign.py",
               "run_targeted_rescue.py", "run_cobalt_carbonyl_search.py", "retry_ionic_reference.py",
               "diagnose_ts_candidate.py"}
    live = []
    for process in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            command = process.info["cmdline"] or []
            if "python" not in (process.info["name"] or "").lower():
                continue
            scripts = [arg for arg in command if str(arg).replace("\\", "/").split("/")[-1] in tracked]
            if scripts and any("pincer-catmech-ai" in arg for arg in command):
                live.append({"pid": process.pid, "scripts": scripts})
        except (psutil.Error, OSError):
            continue
    if live:
        raise RuntimeError(f"Scientific or telemetry writers are still active: {live}")
    neb_summary_path = REPO / "data/campaign/neb/FINAL_NEB_SUMMARY.json"
    neb_summary = load(neb_summary_path)
    if not neb_summary.get("all_tracked_actual_neb_calculation_sessions_ended"):
        raise ValueError("Final NEB identity and evidence audit is incomplete")
    audited_routes = {row["directory"]: row for group in neb_summary["groups"] for row in group["routes"]}
    tasks = []
    for path in sorted((REPO / "data/campaign/neb").glob("*/task.json")):
        data = load(path)
        result_path = path.parent / "path/result.json"
        result = load(result_path) if result_path.exists() else data.get("path_result", {})
        audited = audited_routes.get(path.parent.name)
        if audited is None or digest(path) != audited["task_or_diagnostic_sha256"]:
            raise ValueError(f"NEB task is missing from, or changed after, the final identity audit: {path}")
        if result_path.exists() and digest(result_path) != audited["result_sha256"]:
            raise ValueError(f"NEB result changed after the final identity audit: {result_path}")
        tasks.append({"catalyst_id": data.get("catalyst_id"), "task": path.relative_to(REPO).as_posix(),
                      "task_sha256": digest(path), "status": data.get("status", result.get("status", "unclassified")),
                      "mode": data.get("mode", "search"), "accepted": bool(audited["accepted"]),
                      "engine_accepted": audited["engine_accepted"],
                      "intended_identity_eligible": audited["intended_identity_eligible"],
                      "identity_exclusion": audited["identity_exclusion"],
                      "neb_converged": result.get("neb_converged"),
                      "calculator_evaluations": result.get("calculator_evaluations", 0),
                      "failure_reasons": audited["failure_reasons"]})
    if len(tasks) != len(audited_routes) or sum(row["accepted"] for row in tasks) != neb_summary["accepted_path_count"]:
        raise ValueError("Frozen NEB inventory disagrees with the final identity audit")
    thermo = verify_thermochemistry_snapshot()
    best, rows = freeze_structures()
    conditions = load(REPO / "data/campaign/conditions.json")
    kinetic = campaign_readiness({}, temperatures=conditions["temperature_grid_K"], require_base_path=True)
    ml = training_readiness([])
    # An electronic TS alone would still not supply every required molecular
    # free energy and six reaction channels. No unconstructed network is run.
    atomic_json(REPO / "data/datasets/kinetics_readiness.json", {**kinetic,
                "executed_grid_points": 0, "computed_TOF_count": 0,
                "reason_scope": "No complete source-verified thermodynamic network snapshot was assembled",
                "conditions": conditions})
    atomic_json(REPO / "data/datasets/surrogate_readiness.json", {**ml,
                "trained": False, "model_file": None, "Pareto_front": None,
                "available_validated_barrier_labels": 0})
    atomic_json(REPO / "data/datasets/neb_search_inventory.json", tasks)
    alpb_records = REPO / "data/campaign_alpb_toluene/records_with_identity.jsonl"
    native_records = [json.loads(line) for line in alpb_records.open(encoding="utf-8") if line.strip()]
    successful = [row for row in native_records if row["native_result"].get("converged")]
    state_counts = Counter(record["state"] for record in best.values())
    status = {"schema_version": 1, "frozen_at_utc": utc_now(),
        "campaign_started_utc": original["started_at_utc"],
        "original_deadline_epoch": original["deadline_epoch"],
        "main_finished_utc": alpb["finished_at_utc"],
        "actual_main_window_seconds": alpb["started_epoch"] + alpb["elapsed_seconds"] - original["started_epoch"],
        "publication_status": "prepared_for_verified_publication",
        "scientific_status": "partial_results_with_unmet_transition_state_kinetic_and_surrogate_targets",
        "designs": 24, "intended_state_combinations": 72,
        "gas_preoptimization_jobs": record_count(REPO / "data/campaign/records.jsonl"),
        "alpb_conformer_jobs": len(native_records), "alpb_native_converged_jobs": len(successful),
        "alpb_identity_preserving_converged_jobs": sum(r.get("chemical_identity", {}).get("retained", False) for r in successful),
        "alpb_assembly_failure_records": record_count(REPO / "data/campaign_alpb_toluene/assembly_failures.jsonl"),
        "best_found_states": dict(state_counts), "best_found_state_total": len(best),
        "missing_active_designs": [r["catalyst_id"] for r in rows if f"{r['catalyst_id']}|active" not in best],
        "certified_thermochemistry": {key: thermo[key] for key in ("certificate_count", "eligible_certificate_count", "selected_catalyst_states", "selected_active_catalysts", "selected_reference_species", "catalyst_temperature_rows", "reference_temperature_rows", "reaction_temperature_rows")},
        "thermodynamic_cycle_checks": {key: thermo[key] for key in ("cycle_temperature_rows", "cycle_closure_tolerance_kcal_mol", "cycle_maximum_absolute_residual_kcal_mol", "cycle_interpretation")},
        "neb_search_records": len(tasks), "neb_distinct_catalysts_attempted": sorted({r["catalyst_id"] for r in tasks if r["catalyst_id"]}),
        "accepted_dehydrogenation_TS_records": sum(r["accepted"] for r in tasks),
        "kinetic_grid_points_executed": 0, "trained_surrogate_models": 0,
        "independent_ionic_contact_pair_minimum": load(REPO / "data/ionic_reference_retry/distance2.2_warm1000/result.json")["accepted"],
        "conditions": conditions,
        "interpretation": ["Designed analogues, not a set of experimentally established catalysts",
             "Best found sampled conformers and certified minima, not proven global minima",
             "GFN2 nominal occupations do not determine ground spin",
             "No experimentally validated activity or reaction-controlling step",
             "Failed or wrong-identity paths cannot provide kinetic barriers",
             "TOF/Pareto and trained surrogate objectives remain unmet"],
        "source_sha256": {p.relative_to(REPO).as_posix(): digest(p) for p in
            (alpb_records, REPO / "data/campaign/control.json", REPO / "data/campaign_alpb_toluene/control.json",
             REPO / "data/datasets/catalyst_descriptors.csv", REPO / "data/datasets/thermochemistry/export_protocol.json",
             neb_summary_path)}}
    atomic_json(REPO / "data/CAMPAIGN_STATUS.json", status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    main()
