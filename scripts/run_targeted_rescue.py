"""Additional real searches of Ph-bipyridine designs with lost P coordination.

This independent search preserves the main campaign's immutable attempts and
uses fresh, unrotated ligand embeddings. It never weakens coordination gates.
"""
from dataclasses import asdict
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np
from ase.io import read, write
from pincer_catmech.generators.combinatorial_pincer import enumerate_specs, generate_catalyst
from pincer_catmech.features.steric import bite_angle, buried_volume
from pincer_catmech.quantum.xtb_backend import optimize_geometry, atomic_json, utc_now
from pincer_catmech.quantum.thermochemistry import characterize_stationary_point
from run_reference_thermochemistry import connectivity_check


def ligand_connectivity(atoms, metadata):
    indices = metadata["ligand_indices"]
    remap = {old: new for new, old in enumerate(indices)}
    bonds = [(remap[i], remap[j]) for i, j in metadata["ligand_bond_indices"]]
    return connectivity_check(atoms[indices], bonds)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deadline-epoch", required=True, type=float)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    data = REPO / "data" / "targeted_rescue"
    scratch = REPO.parent / "targeted-rescue-scratch"
    data.mkdir(parents=True, exist_ok=True)
    records, best = [], {}
    for state in ("active", "hydrogenated", "protonated_reference"):
        for spec in enumerate_specs():
            if spec.backbone != "bipyridine_pnnoh" or spec.substituent != "Ph":
                continue
            for index, seed in enumerate((20260913, 20270000, 20270817, 20273333)):
                if time.time() >= args.deadline_epoch:
                    break
                identifier = f"{spec.catalyst_id}__{state}__seed{seed}"
                folder = scratch / identifier
                try:
                    generated = generate_catalyst(spec, state, seed=seed, arm_torsion_deg=0)
                    metadata = generated.metadata()
                    native = optimize_geometry(generated.atoms, folder, charge=metadata["charge"],
                        unpaired=0, threads=args.threads, max_cycles=400, convergence="vtight",
                        solvent="toluene", timeout_seconds=min(600, max(1, args.deadline_epoch - time.time())))
                    record = {"catalyst_id": spec.catalyst_id, "state": state, "seed": seed,
                        "conformer": f"rescue_seed{seed}", "native_result": asdict(native),
                        "metadata": metadata, "solvent": "toluene", "solvation_state": "gsolv",
                        "accepted_geometry_screen": False, "evaluated_at_utc": utc_now()}
                    if native.converged:
                        atoms = read(folder / "xtbopt.xyz")
                        distances = atoms.get_distances(0, metadata["donor_indices"])
                        connectivity = ligand_connectivity(atoms, metadata)
                        retained = bool(np.all((distances > 1.4) & (distances < 3.1)))
                        record.update(metal_donor_distances_A=distances.tolist(),
                            coordination_retained=retained, ligand_connectivity=connectivity,
                            accepted_geometry_screen=retained and connectivity["retained"],
                            optimized_geometry_sha256=hashlib.sha256((folder / "xtbopt.xyz").read_bytes()).hexdigest())
                        if record["accepted_geometry_screen"]:
                            volume = buried_volume(atoms, 0, ligand_indices=metadata["ligand_indices"], n_samples=50000, seed=2026)
                            record.update(buried_volume_percent=volume.percent_buried_volume,
                                buried_volume_Wilson95_percent=list(volume.confidence_interval_percent),
                                terminal_donor_angle_deg=bite_angle(atoms, 0, metadata["donor_indices"][0], metadata["donor_indices"][2]),
                                donor_labels="P-N-N", P_M_P_angle_deg=None)
                            key = f"{spec.catalyst_id}|{state}"
                            if key not in best or native.energy_eV < best[key]["native_result"]["energy_eV"]:
                                best[key] = record
                    atomic_json(folder / "record.json", record)
                except Exception as exc:
                    record = {"catalyst_id": spec.catalyst_id, "state": state, "seed": seed,
                              "status": "failed", "error": f"{type(exc).__name__}: {exc}", "utc": utc_now()}
                records.append(record)
                with (data / "records.jsonl").open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(record, allow_nan=False) + "\n")
                atomic_json(data / "best_conformers.json", best)
                print(json.dumps({"utc": utc_now(), "identifier": identifier,
                    "accepted_geometry": record.get("accepted_geometry_screen", False),
                    "distances_A": record.get("metal_donor_distances_A"), "error": record.get("error")}), flush=True)
    for key, record in best.items():
        if time.time() >= args.deadline_epoch:
            break
        name = key.replace("|", "__")
        atoms = read(Path(record["native_result"]["output_directory"]) / "xtbopt.xyz")
        result = characterize_stationary_point(atoms, data / "thermochemistry" / name,
            scratch / "thermochemistry" / name, charge=record["metadata"]["charge"],
            threads=args.threads, max_repairs=2, deadline_epoch=args.deadline_epoch,
            coordination_indices=record["metadata"], source_metadata=record, solvent="toluene")
        print(json.dumps({"utc": utc_now(), "stationarity": name, "status": result["status"]}), flush=True)
    atomic_json(data / "summary.json", {"finished_utc": utc_now(), "attempts": len(records),
                "accepted_species": len(best), "method": "GFN2-xTB/ALPB toluene gsolv",
                "scope": "Fresh unrotated embedding seeds; unrestrained native optimizations; additional local search"})


if __name__ == "__main__":
    main()
