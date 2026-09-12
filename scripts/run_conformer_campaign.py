"""Continuous, checkpointed native GFN2-xTB conformer exploration.

Each submitted geometry is independently embedded/torsion-perturbed. Search
continues until the requested wall-clock deadline; elapsed time is measured,
not simulated. Every completed or failed native run has its own provenance.
"""

from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import psutil
from ase.io import read, write

from pincer_catmech.features.steric import bite_angle, buried_volume
from pincer_catmech.generators.combinatorial_pincer import enumerate_specs, generate_catalyst
from pincer_catmech.quantum.xtb_backend import optimize_geometry, OptimizationResult, atomic_json, utc_now
from pincer_catmech.quantum.identity import catalyst_identity_screen


def campaign_protocol(solvent):
    """The fixed scientific protocol; resource settings are deliberately separate.

    gsolv is the configured xTB solvation reference convention. For gas-phase
    jobs no solvation argument is passed to xTB, as recorded by solvent.
    """
    if solvent not in (None, "toluene"):
        raise ValueError("The campaign supports gas phase or ALPB toluene")
    return {"method": "GFN2-xTB", "electronic_temperature_K": 300.0,
            "solvent": solvent or "gas phase", "solvation_state": "gsolv"}


def validate_resume_protocol(previous, solvent):
    """Reject mixed potentials before modifying control, records, or selections."""
    expected = campaign_protocol(solvent)
    mismatches = [f"{field}: recorded={previous.get(field)!r}, requested={value!r}"
                  for field, value in expected.items() if previous.get(field) != value]
    if mismatches:
        raise ValueError("Cannot resume with a different or incomplete scientific protocol; "
                         "start a new campaign name. " + "; ".join(mismatches))


def optimize_job(spec, state, conformer, data_root, executable, threads, deadline, solvent=None):
    root = Path(data_root)
    seed = 17000 + conformer
    generated = generate_catalyst(spec, state=state, seed=seed,
                                 arm_torsion_deg=((conformer * 37) % 121) - 60)
    metadata = generated.metadata()
    identifier = f"{spec.metal}_{spec.backbone}_{spec.substituent}"
    run_family = "campaign-runs" if solvent is None else f"campaign-runs-alpb-{solvent}"
    if root.name not in ("campaign", "campaign_alpb_toluene"):
        run_family = "campaign-runs-" + root.name
    folder = REPO.parent / run_family / identifier / state / f"conformer_{conformer:04d}"
    folder.mkdir(parents=True, exist_ok=True)
    atomic_json(folder / "structure_metadata.json", metadata)
    write(folder / "initial.xyz", generated.atoms)
    if (folder / "result.json").exists():
        result = OptimizationResult(**json.loads((folder / "result.json").read_text()))
    else:
        result = optimize_geometry(generated.atoms, folder, charge=metadata["charge"],
            unpaired=metadata["unpaired_electrons"], threads=threads, executable=executable,
            max_cycles=300, timeout_seconds=max(1, min(420, deadline-time.time())), solvent=solvent)
    record = {"catalyst_id": identifier, "metal": spec.metal, "backbone": spec.backbone,
        "substituent": spec.substituent, "state": state, "conformer": conformer,
        "seed": seed, "native_result": asdict(result), "metadata": metadata,
        "evaluated_at_utc": utc_now(), "coordination_retained": False,
        "solvent": solvent, "solvation_state": "gsolv"}
    if result.converged:
        atoms = read(folder / "xtbopt.xyz")
        donors = tuple(metadata["donor_indices"])
        distances = atoms.get_distances(0, list(donors))
        retained = bool(np.all((distances > 1.4) & (distances < 3.1)))
        record["coordination_retained"] = retained
        record["chemical_identity"] = catalyst_identity_screen(atoms, metadata)
        record["metal_donor_distances_A"] = distances.tolist()
        record["terminal_donor_angle_deg"] = bite_angle(atoms, 0, donors[0], donors[2])
        record["donor_labels"] = "-".join(atoms[i].symbol for i in donors)
        record["P_M_P_angle_deg"] = (record["terminal_donor_angle_deg"]
                                    if atoms[donors[0]].symbol == atoms[donors[2]].symbol == "P" else None)
        volume = buried_volume(atoms, 0, ligand_indices=metadata["ligand_indices"],
                               n_samples=50000, seed=2026)
        record["buried_volume_percent"] = volume.percent_buried_volume
        record["buried_volume_Wilson95_percent"] = list(volume.confidence_interval_percent)
        record["optimized_geometry_sha256"] = hashlib.sha256((folder / "xtbopt.xyz").read_bytes()).hexdigest()
    atomic_json(folder / "record.json", record)
    return record


def write_descriptors(data_root, records):
    best, tried = {}, {}
    for record in records:
        key = (record["catalyst_id"], record["state"])
        tried[key] = tried.get(key, 0) + 1
        if (record["native_result"]["converged"] and record["coordination_retained"]
                and record.get("chemical_identity", {}).get("retained") is True):
            if key not in best or record["native_result"]["energy_eV"] < best[key]["native_result"]["energy_eV"]:
                best[key] = record
    fields = ["catalyst_id", "metal", "backbone", "substituent", "state", "status",
              "conformers_attempted", "best_conformer", "energy_eV", "max_force_eV_A",
              "terminal_donor_angle_deg", "donor_labels", "P_M_P_angle_deg", "buried_volume_percent",
              "buried_volume_CI_low_percent", "buried_volume_CI_high_percent", "geometry_path",
              "proton_site_element", "charge", "unpaired_electrons"]
    rows = []
    for spec in enumerate_specs():
        identifier = f"{spec.metal}_{spec.backbone}_{spec.substituent}"
        key = (identifier, "active")
        record = best.get(key)
        row = {"catalyst_id": identifier, "metal": spec.metal, "backbone": spec.backbone,
               "substituent": spec.substituent, "state": "active", "conformers_attempted": tried.get(key, 0),
               "status": "best_found_identity_preserved_optimized_geometry" if record else "no_accepted_conformer"}
        if record:
            row.update({name: record[name] for name in ("terminal_donor_angle_deg", "donor_labels",
                                                       "P_M_P_angle_deg", "buried_volume_percent")})
            row.update({"best_conformer": record["conformer"],
                "energy_eV": record["native_result"]["energy_eV"],
                "max_force_eV_A": record["native_result"]["max_force_eV_A"],
                "buried_volume_CI_low_percent": record["buried_volume_Wilson95_percent"][0],
                "buried_volume_CI_high_percent": record["buried_volume_Wilson95_percent"][1],
                "geometry_path": str(Path(record["native_result"]["output_directory"]) / "xtbopt.xyz"),
                "proton_site_element": record["metadata"]["proton_site_element"],
                "charge": record["metadata"]["charge"], "unpaired_electrons": record["metadata"]["unpaired_electrons"]})
        rows.append(row)
    output = (REPO / "data" / "datasets" if Path(data_root).name in ("campaign", "campaign_alpb_toluene")
              else Path(data_root) / "datasets")
    output.mkdir(parents=True, exist_ok=True)
    with (output / "catalyst_descriptors.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    atomic_json(output / "best_conformers.json", {"|".join(key): val for key, val in best.items()})
    return best


def run_campaign(hours, workers, threads, executable, resume=False, solvent=None, deadline_epoch=None, campaign_name=None):
    if not 0 < hours <= 24 or not 1 <= workers <= 8 or not 1 <= threads <= 8:
        raise ValueError("Invalid campaign resource bounds")
    if campaign_name is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", campaign_name):
        raise ValueError("campaign name must contain only letters, digits, underscore and hyphen")
    protocol = campaign_protocol(solvent)
    data_root = REPO / "data" / (campaign_name or ("campaign" if solvent is None else f"campaign_alpb_{solvent}"))
    data_root.mkdir(parents=True, exist_ok=True)
    control = data_root / "control.json"
    previous = json.loads(control.read_text()) if control.exists() else None
    if resume and previous is None:
        raise ValueError("Cannot resume without an existing control.json and its original scientific protocol")
    if previous is not None and not isinstance(previous, dict):
        raise ValueError("Existing campaign control must be a protocol-bearing JSON object")
    if previous is not None and previous.get("status") == "running" and psutil.pid_exists(previous["pid"]):
        raise RuntimeError("A campaign process is already running")
    if previous is not None and not resume:
        raise RuntimeError("Campaign records exist; use --resume to preserve its original clock and calculations")
    if previous is not None:
        validate_resume_protocol(previous, solvent)
    start = previous["started_epoch"] if previous else time.time()
    deadline = previous["deadline_epoch"] if previous else (deadline_epoch or start + hours * 3600)
    info = {"status": "running", "started_at_utc": previous["started_at_utc"] if previous else utc_now(), "started_epoch": start,
            "deadline_epoch": deadline, "requested_hours": previous["requested_hours"] if previous else hours,
            "pid": os.getpid(), "workers": workers, "threads_per_worker": threads, **protocol,
            "standard_state_correction": "added separately in qRRHO at each temperature",
            "interpretation": "Designed analogues; best found structures, no global-minimum or ground-spin proof"}
    atomic_json(control, info)
    records = [json.loads(line) for line in (data_root / "records.jsonl").read_text().splitlines()] if (data_root / "records.jsonl").exists() else []
    # Earlier raw records are immutable. A derived identity-audited view fixes
    # the initially insufficient donor-distance-only selection criterion.
    for record in records:
        if record["native_result"]["converged"]:
            geometry = Path(record["native_result"]["output_directory"]) / "xtbopt.xyz"
            if hashlib.sha256(geometry.read_bytes()).hexdigest() != record["optimized_geometry_sha256"]:
                raise ValueError("Stored optimized-geometry checksum changed")
            record["chemical_identity"] = catalyst_identity_screen(read(geometry), record["metadata"])
    with (data_root / "records_with_identity.jsonl").open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
    if records:
        write_descriptors(data_root, records)
    seen = {(r["catalyst_id"], r["state"], r["conformer"]) for r in records}
    pending = {}
    specs = enumerate_specs()
    round_number, spec_index, state_index = 0, 0, 0
    states = ("active", "hydrogenated", "protonated_reference")
    last_telemetry = start - 600
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        while time.time() < deadline or pending:
            while len(pending) < workers and time.time() < deadline:
                spec, state = specs[spec_index], states[state_index]
                conformer = round_number
                spec_index += 1
                if spec_index == len(specs):
                    spec_index = 0
                    state_index += 1
                    if state_index == len(states):
                        state_index = 0
                        round_number += 1
                if (spec.catalyst_id, state, conformer) in seen:
                    continue
                future = pool.submit(optimize_job, spec, state, conformer,
                                     str(data_root), executable, threads, deadline, solvent)
                pending[future] = (spec.metal, spec.backbone, spec.substituent, state, conformer)
            completed, _ = wait(pending, timeout=2, return_when=FIRST_COMPLETED)
            for future in completed:
                label = pending.pop(future)
                try:
                    record = future.result()
                    key = (record["catalyst_id"], record["state"], record["conformer"])
                    if key not in seen:
                        records.append(record)
                        seen.add(key)
                        with (data_root / "records.jsonl").open("a", encoding="utf-8") as stream:
                            stream.write(json.dumps(record, allow_nan=False) + "\n")
                        with (data_root / "records_with_identity.jsonl").open("a", encoding="utf-8") as stream:
                            stream.write(json.dumps(record, allow_nan=False) + "\n")
                except Exception as exc:
                    failures += 1
                    with (data_root / "assembly_failures.jsonl").open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps({"label": label, "utc": utc_now(), "error": repr(exc)}) + "\n")
            if completed:
                best = write_descriptors(data_root, records)
                info.update(completed_jobs=len(records), assembly_failures=failures,
                            accepted_species=len(best), active_jobs=len(pending), current_round=round_number,
                            elapsed_seconds=time.time()-start)
                atomic_json(control, info)
            if time.time() - last_telemetry >= 600:
                last_telemetry = time.time()
                processes = []
                for child in psutil.Process().children(recursive=True):
                    try:
                        processes.append({"pid": child.pid, "name": child.name(),
                            "cpu_seconds": sum(child.cpu_times()[:2]), "rss_bytes": child.memory_info().rss})
                    except psutil.Error:
                        continue
                telemetry = {"utc": utc_now(), "elapsed_seconds": time.time()-start,
                    "completed_jobs": len(records), "native_converged_jobs": sum(r["native_result"]["converged"] for r in records),
                    "assembly_failures": failures, "active_jobs": list(pending.values()),
                    "local_compute_processes": processes, "host_cpu_percent": psutil.cpu_percent(),
                    "last_energy_eV": records[-1]["native_result"]["energy_eV"] if records else None,
                    "hosted_agent_cpu": "not observable; CPU metrics above belong to local scientific workers"}
                with (data_root / "telemetry.jsonl").open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(telemetry) + "\n")
                print(json.dumps(telemetry), flush=True)
    info.update(status="finished", finished_at_utc=utc_now(), elapsed_seconds=time.time()-start,
                completed_jobs=len(records), assembly_failures=failures)
    atomic_json(control, info)
    write_descriptors(data_root, records)
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hours", type=float, default=3)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--xtb", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--solvent", choices=["toluene"])
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--campaign-name", help="Fresh run identifier, preserving published data and using a separate scratch directory")
    options = parser.parse_args()
    run_campaign(options.hours, options.workers, options.threads, options.xtb,
                 options.resume, options.solvent, options.deadline_epoch, options.campaign_name)
