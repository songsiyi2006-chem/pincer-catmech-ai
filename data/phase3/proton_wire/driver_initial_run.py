"""Bounded physical neutral tBuOH-assisted dehydration CI-NEB exploration.

Uses mapped Phase3 open-cluster minima and bond-order-correct accepted imine,
water, and tBuOH component geometries. No potassium/base-order claim is made.
All native evidence, including failed endpoints and incomplete bands, is kept.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import time
import zipfile

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


def archive_evidence(scratch, destination):
    entries = []
    with zipfile.ZipFile(destination / "native_evidence.zip", "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(Path(scratch).rglob("*")):
            if path.is_file():
                member = path.relative_to(scratch).as_posix()
                archive.write(path, member)
                entries.append({"member": member, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
    artifact = provenance(destination / "native_evidence.zip")
    atomic_json(destination / "native_manifest.json", {"archive": artifact, "members": entries,
        "member_count": len(entries), "raw_bytes": sum(row["bytes"] for row in entries),
        "scope": "Every native file within this new campaign scratch tree; prior campaigns unchanged"})
    return {"archive": artifact, "member_count": len(entries)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb", required=True)
    parser.add_argument("--budget-seconds", type=float, default=1200)
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--neb-steps", type=int, default=100)
    args = parser.parse_args()
    start = time.time()
    if not 60 <= args.budget_seconds <= 1200 or args.neb_steps < 1:
        raise ValueError("Require a bounded 60..1200 second budget and positive NEB steps")
    deadline = min(start + args.budget_seconds, args.deadline_epoch or float("inf"))
    output = REPO / "data" / "phase3" / "proton_wire"
    scratch = REPO.parent / "phase3" / "proton_wire"
    if output.exists() or scratch.exists():
        raise FileExistsError("Fresh output and scratch paths required; preserve completed evidence")
    output.mkdir(parents=True)
    scratch.mkdir(parents=True)
    summary = {"started_utc": utc_now(), "status": "running", "method": "GFN2-xTB", "solvent": "toluene",
        "solvation_state": "gsolv", "physical_workers": 1, "threads": 1, "charge": 0, "unpaired_electrons": 0,
        "thermodynamic_baseline_temperature_K": 383.15, "electronic_smearing_temperature_K": 300.0,
        "tBuOK_total_base_equivalent_baseline": .05, "total_base_speciation_resolved": False,
        "potassium_explicit": False, "tert_butoxide_explicit": False,
        "scope": "Neutral micro-solvated dehydration hypothesis; no base kinetic order or activation free energy inferred",
        "deadline_epoch": deadline, "attempts": [], "accepted_TS_count": 0, "activation_free_energy_barrier_count": 0}
    atomic_json(output / "summary.json", summary)
    try:
        for number, (seed, distance) in enumerate(((2, 4.3), (1, 5.4))):
            if deadline-time.time() < 90:
                summary["unattempted_second_reason"] = "Remaining bounded wall time below 90 seconds"
                break
            tag = f"attempt_{number+1:02d}_cluster_seed_{seed:02d}"
            prepared = output / (tag + "_assembly")
            prepared.mkdir()
            source = REPO / "data" / "phase3" / "solvation" / "clusters" / f"hemiaminal_tbuoh_1_seed_{seed:02d}"
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
        summary["native_evidence"] = archive_evidence(scratch, output)
        atomic_json(output / "summary.json", summary)
        print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
