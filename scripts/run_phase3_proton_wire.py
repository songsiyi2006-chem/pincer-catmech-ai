"""Bounded physical neutral tBuOH-assisted dehydration CI-NEB exploration.

Uses mapped Phase3 open-cluster minima and bond-order-correct accepted imine,
water, and tBuOH component geometries. No potassium/base-order claim is made.
All native evidence, including failed endpoints and incomplete bands, is kept.
Use --finalize-evidence after the initial run and optional continuation to make
verified complete ZIP parts below 95 MiB. Finalization never runs quantum jobs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from datetime import datetime

for _thread_variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np
from ase.io import read, write
from pincer_catmech.kinetics.condensation_search import identity, observed_bonds, run_condensation_path
from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now
from run_condensation_search import product_from_reference_components
from finalize_phase3_proton_wire import finalize_evidence, separate_paths, verify_archives


def provenance(path):
    path = Path(path)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size}


def mapped_endpoint_seed(cluster_dir, water_distance=4.3):
    """Prepare chemically annotated seeds without evaluating an energy."""
    cluster_dir = Path(cluster_dir)
    metadata_path = cluster_dir / "mapping.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata["charge"] != 0 or len(metadata["fragments"]) != 2:
        raise ValueError("This protocol requires a neutral substrate plus one tBuOH")
    stationarity = cluster_dir / "stationarity" / "result.json"
    accepted = False
    if stationarity.exists():
        record = json.loads(stationarity.read_text(encoding="utf-8"))
        accepted = record.get("accepted") is True and record.get("status") == "verified_minimum" and record.get("solvent") == "toluene"
    source = cluster_dir / "stationarity" / "accepted_geometry.xyz" if accepted else cluster_dir / "optimized.xyz"
    reactant = read(source)
    bonds = {tuple(sorted(pair)) for pair in metadata["bonds"]}
    atom_map = metadata["atom_map"]
    if [row["global_index"] for row in atom_map] != list(range(len(reactant))):
        raise ValueError("Phase3 source lacks complete ordered atom mapping")
    if [row["symbol"] for row in atom_map] != reactant.get_chemical_symbols():
        raise ValueError("Mapped atomic identities differ from source geometry")
    if observed_bonds(reactant) != bonds:
        raise ValueError("Phase3 source connectivity differs from the annotated reactant")
    sites = metadata["substrate_sites"]
    oxygen, ho, nitrogen, hn = [sites[key] for key in ("oxygen", "hydroxyl_hydrogen", "nitrogen", "amine_hydrogen")]
    carbon = 1
    fragment = metadata["fragments"][1]
    ot = next(i for i in fragment if reactant[i].symbol == "O")
    ht = next(i for i in fragment if reactant[i].symbol == "H" and tuple(sorted((ot, i))) in bonds)
    product = reactant.copy()
    direction = reactant.positions[oxygen] - reactant.positions[carbon]
    direction /= np.linalg.norm(direction)
    product.positions[[oxygen, ho]] += (water_distance - reactant.get_distance(oxygen, carbon)) * direction
    direction = reactant.positions[nitrogen] - reactant.positions[ot]
    product.positions[hn] = product.positions[ot] + direction / np.linalg.norm(direction)
    direction = product.positions[ot] - product.positions[oxygen]
    product.positions[ht] = product.positions[oxygen] + .98 * direction / np.linalg.norm(direction)
    product, references = product_from_reference_components(product, anionic=False, water_distance=water_distance)
    product_bonds = set(bonds)
    broken = [(oxygen, carbon), (nitrogen, hn), (ot, ht)]
    formed = [(ot, hn), (oxygen, ht)]
    for pair in broken:
        product_bonds.remove(tuple(sorted(pair)))
    for pair in formed:
        product_bonds.add(tuple(sorted(pair)))
    coordinates = {
        "N_H_to_tBuOH_O": {"terms": [(nitrogen, hn, 1), (ot, hn, -1)], "minimum_overlap": .10},
        "tBuOH_H_to_leaving_O": {"terms": [(ot, ht, 1), (oxygen, ht, -1)], "minimum_overlap": .10},
        "C_O_cleavage": {"terms": [(carbon, oxygen, 1)], "minimum_overlap": .05},
        "C_N_contraction": {"terms": [(carbon, nitrogen, -1)], "minimum_overlap": .03},
    }
    details = {"source_geometry": provenance(source), "source_mapping": provenance(metadata_path),
        "source_already_hessian_certified": accepted, "source_stationarity": provenance(stationarity) if stationarity.exists() else None,
        "product_component_references": references, "atom_map": atom_map,
        "mapped_indices": {"C": carbon, "N": nitrogen, "leaving_O": oxygen, "original_O_H": ho,
                           "N_H": hn, "shuttle_O": ot, "shuttle_H": ht},
        "broken_bonds": broken, "formed_bonds": formed, "C_N_bond_order_change": "single_to_double",
        "product_seed_C_O_distance_A": water_distance,
        "product_seed_identity": identity(product, product_bonds, cn=(carbon, nitrogen), cn_interval=(1.15, 1.42)),
        "reactant_seed_identity": identity(reactant, bonds, cn=(carbon, nitrogen), cn_interval=(1.35, 1.70)),
        "seed_geometry_is_stationary_point_evidence": False,
        "proton_wire": "Neutral tBuOH receives original N-H while its original O-H protonates leaving water oxygen",
        "closed_six_or_eight_membered_TS_demonstrated": False,
        "potassium_explicit": False, "tert_butoxide_explicit": False}
    return reactant, product, [sorted(bonds), sorted(product_bonds)], coordinates, (carbon, nitrogen), details


def report_attempt(result, directory):
    snapshots = [stage for stage in result.get("stages", []) if stage["stage"] == "neb"]
    energies = result.get("image_energies_eV") or (snapshots[-1]["image_energies_eV"] if snapshots else None)
    accepted = result.get("accepted") is True and result.get("status") == "verified_first_order_saddle"
    return {"directory": str(directory), "status": result["status"], "accepted_TS": accepted,
        "failure_reasons": result.get("failure_reasons", []),
        "actual_NEB_evaluation_reached": bool(snapshots), "last_reported_NEB_step": snapshots[-1]["step"] if snapshots else None,
        "neb_converged": result.get("neb_converged", False),
        "last_sampled_band_potential_max_above_reactant_eV": float(max(energies)-energies[0]) if energies else None,
        "sampled_band_potential_is_activation_free_energy": False,
        "accepted_endpoint_minima": sum(row.get("stationarity_status") == "verified_minimum" for row in result.get("endpoints", [])),
        "imaginary_frequency_count": result.get("imaginary_count"),
        "free_energy_barrier_computed": result.get("free_energy_barrier_computed", False),
        "elapsed_seconds": result.get("elapsed_seconds"), "result": provenance(directory / "result.json")}


def execution_source(output, stage):
    """Snapshot this exact prospective driver before a physical stage starts."""
    destination = output / "run_sources" / (stage + "_driver.py")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError("Executed source snapshots are immutable")
    destination.write_bytes(Path(__file__).read_bytes())
    return {"stage": stage, "relative_path": destination.relative_to(output).as_posix(),
        "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(), "captured_utc": utc_now()}


def pending_evidence(scratch):
    return {"status": "pending_finalization", "raw_files_preserved": True,
        "scratch_directory": str(scratch), "next_action": "Run --finalize-evidence after the optional continuation"}


def continue_first_attempt(args):
    """One immutable continuation under the original campaign wall deadline."""
    output, scratch = separate_paths(args.output, args.scratch)
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    if (output / "native_manifest.json").exists() or (output / "verification.json").exists():
        raise FileExistsError("Finalized evidence is immutable; continue only before --finalize-evidence")
    if summary["status"] == "running" or summary.get("continuation_count", 0):
        raise ValueError("Require completed initial attempts and no existing continuation")
    deadline = min(summary["deadline_epoch"], args.deadline_epoch or float("inf"))
    if deadline-time.time() < 60:
        raise TimeoutError("Original physical campaign deadline leaves under 60 seconds")
    prior = output / "attempt_01_cluster_seed_02"
    previous = json.loads((prior / "result.json").read_text(encoding="utf-8"))
    protocol = json.loads((prior / "protocol.json").read_text(encoding="utf-8"))
    if previous["status"] == "running" or previous.get("accepted"):
        raise ValueError("Continuation requires a completed unaccepted band")
    band = prior / "final_band.xyz"
    if not band.exists():
        band = prior / "current_band.xyz"
    frames = read(band, index=":")
    if len(frames) != 9 or any(not np.array_equal(frame.numbers, frames[0].numbers) for frame in frames):
        raise ValueError("Restart must contain nine identically mapped images")
    tag = "attempt_01_continuation"
    prepared = output / (tag + "_assembly")
    prepared.mkdir(exist_ok=False)
    metadata = json.loads((output / (prior.name + "_assembly") / "mapping_and_provenance.json").read_text(encoding="utf-8"))
    metadata.update(restart_source=provenance(band), restart_protocol=provenance(prior / "protocol.json"),
        previous_result=provenance(prior / "result.json"), new_interpolation=False, all_forces_recomputed=True)
    atomic_json(prepared / "mapping_and_provenance.json", metadata)
    history = output / "before_continuation_summary.json"
    if history.exists():
        raise FileExistsError("Previous summary history is immutable")
    shutil.copy2(output / "summary.json", history)
    summary.setdefault("executed_driver_sources", []).append(execution_source(output, "continuation"))
    summary.update(status="running", continuation_count=1, distinct_seed_attempt_count=2,
        continuation_started_utc=utc_now(), original_physical_deadline_preserved=True)
    atomic_json(output / "summary.json", summary)
    try:
        result = run_condensation_path(frames[0], frames[-1], output_dir=output / tag, scratch_dir=scratch / tag,
            executable=args.xtb, charge=0, endpoint_bonds=protocol["endpoint_bonds"],
            reaction_coordinates=protocol["reaction_coordinates"], cn=protocol["C_N_indices"],
            deadline_epoch=deadline, neb_steps=args.neb_steps, endpoint_cycles=500,
            neb_fmax=protocol["neb_fmax_eV_A"], ts_fmax=protocol["ts_fmax_eV_A"],
            neb_maxstep=protocol["neb_maxstep_angstrom"], certify_endpoints=True, restart_band=band,
            model_label=protocol["model"] + "; immutable continuation with fresh forces and no interpolation reset")
        summary["attempts"].append(report_attempt(result, output / tag))
        summary["accepted_TS_count"] = sum(row["accepted_TS"] for row in summary["attempts"])
        summary["activation_free_energy_barrier_count"] = sum(row["free_energy_barrier_computed"] for row in summary["attempts"])
        summary["status"] = "completed_bounded_exploration"
    except Exception as error:
        summary.update(status="continuation_driver_failed", driver_failure=f"{type(error).__name__}: {error}")
    finally:
        started_epoch = summary.get("started_epoch", datetime.fromisoformat(summary["started_utc"]).timestamp())
        summary.update(finished_utc=utc_now(), elapsed_seconds=time.time()-started_epoch)
        summary["native_evidence"] = pending_evidence(scratch)
        atomic_json(output / "summary.json", summary)
        print(json.dumps(summary, ensure_ascii=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb", help="Required only for a new physical run or continuation")
    parser.add_argument("--output", type=Path, default=REPO / "data" / "phase3" / "proton_wire")
    parser.add_argument("--scratch", type=Path, default=REPO.parent / "phase3" / "proton_wire")
    parser.add_argument("--cluster-source", type=Path, default=REPO / "data" / "phase3" / "solvation" / "clusters",
        help="Phase3 mapped cluster directory; allows a fresh master-campaign solvation output")
    parser.add_argument("--budget-seconds", type=float, default=1200)
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--neb-steps", type=int, default=100)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--continue-first", action="store_true")
    mode.add_argument("--finalize-evidence", action="store_true", help="Archive/verify/tabulate completed runs; no QC calculations")
    mode.add_argument("--verify-evidence", action="store_true", help="Verify published ZIPs read-only; no receipt changes")
    args = parser.parse_args()
    output, scratch = separate_paths(args.output, args.scratch)
    if args.finalize_evidence:
        print(json.dumps(finalize_evidence(output, scratch), indent=2), flush=True)
        return
    if args.verify_evidence:
        print(json.dumps(verify_archives(output), indent=2), flush=True)
        return
    if args.xtb is None:
        parser.error("--xtb is required for physical calculations")
    if args.neb_steps < 1 or (args.deadline_epoch is not None and not np.isfinite(args.deadline_epoch)):
        parser.error("Require positive NEB steps and a finite explicit deadline")
    if args.continue_first:
        continue_first_attempt(args)
        return
    start = time.time()
    if not 60 <= args.budget_seconds <= 1200 or args.neb_steps < 1:
        raise ValueError("Require a bounded 60..1200 second budget and positive NEB steps")
    deadline = min(start + args.budget_seconds, args.deadline_epoch or float("inf"))
    if output.exists() or scratch.exists():
        raise FileExistsError("Fresh output and scratch paths required; preserve completed evidence")
    output.mkdir(parents=True)
    scratch.mkdir(parents=True)
    summary = {"started_utc": utc_now(), "started_epoch": start, "status": "running", "method": "GFN2-xTB", "solvent": "toluene",
        "solvation_state": "gsolv", "physical_workers": 1, "threads": 1, "charge": 0, "unpaired_electrons": 0,
        "thermodynamic_baseline_temperature_K": 383.15, "electronic_smearing_temperature_K": 300.0,
        "tBuOK_total_base_equivalent_baseline": .05, "total_base_speciation_resolved": False,
        "potassium_explicit": False, "tert_butoxide_explicit": False,
        "scope": "Neutral micro-solvated dehydration hypothesis; no base kinetic order or activation free energy inferred",
        "deadline_epoch": deadline, "attempts": [], "accepted_TS_count": 0, "activation_free_energy_barrier_count": 0,
        "executed_driver_sources": [execution_source(output, "initial")]}
    atomic_json(output / "summary.json", summary)
    try:
        for number, (seed, distance) in enumerate(((2, 4.3), (1, 5.4))):
            if deadline-time.time() < 90:
                summary["unattempted_second_reason"] = "Remaining bounded wall time below 90 seconds"
                break
            tag = f"attempt_{number+1:02d}_cluster_seed_{seed:02d}"
            prepared = output / (tag + "_assembly")
            prepared.mkdir()
            source = args.cluster_source.resolve() / f"hemiaminal_tbuoh_1_seed_{seed:02d}"
            reactant, product, bonds, coordinates, cn, metadata = mapped_endpoint_seed(source, distance)
            atomic_json(prepared / "mapping_and_provenance.json", metadata)
            write(prepared / "reactant_seed.xyz", reactant)
            write(prepared / "product_seed.xyz", product)
            attempt_deadline = min(deadline, time.time()+min(650, args.budget_seconds*.55)) if number == 0 else deadline
            result = run_condensation_path(reactant, product, output_dir=output / tag, scratch_dir=scratch / tag,
                executable=args.xtb, charge=0, endpoint_bonds=bonds, reaction_coordinates=coordinates, cn=cn,
                deadline_epoch=attempt_deadline, neb_steps=args.neb_steps, endpoint_cycles=500,
                neb_fmax=.07, ts_fmax=.03, neb_maxstep=.06, certify_endpoints=True,
                model_label="Phase3 mapped 43-atom neutral hemiaminal/tBuOH dehydration; product references fitted; complete free endpoint and saddle gates")
            summary["attempts"].append(report_attempt(result, output / tag))
            summary["accepted_TS_count"] = sum(row["accepted_TS"] for row in summary["attempts"])
            summary["activation_free_energy_barrier_count"] = sum(row["free_energy_barrier_computed"] for row in summary["attempts"])
            atomic_json(output / "summary.json", summary)
            if result.get("accepted") is True:
                break
        summary["status"] = "completed_bounded_exploration"
    except Exception as error:
        summary.update(status="driver_failed", driver_failure=f"{type(error).__name__}: {error}")
    finally:
        summary.update(finished_utc=utc_now(), elapsed_seconds=time.time()-start)
        summary["native_evidence"] = pending_evidence(scratch)
        atomic_json(output / "summary.json", summary)
        print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
