"""Test alternate Co carbonyl placements after observed ligand-O/CO coupling.

All placements are only starting geometries. Every accepted structure must
survive unrestricted native optimization and the unchanged full identity gate.
"""
from pathlib import Path
from dataclasses import asdict
import argparse
import hashlib
import json
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import numpy as np
from ase.io import read
from pincer_catmech.generators.combinatorial_pincer import CatalystSpec, generate_catalyst
from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.xtb_backend import optimize_geometry, atomic_json, utc_now
from pincer_catmech.quantum.thermochemistry import characterize_stationary_point
from pincer_catmech.features.steric import bite_angle, buried_volume


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deadline-epoch", required=True, type=float)
    args = parser.parse_args()
    data = REPO / "data" / "cobalt_carbonyl_rescue"
    scratch = REPO.parent / "cobalt-carbonyl-scratch"
    data.mkdir(parents=True, exist_ok=True)
    best = {}
    attempts = 0
    vectors = [(0, 1, 0), (0, -1, 0), (0, 0, 1), (1, 1, 1), (-1, 1, 1), (0, 1, -1)]
    for substituent in ("Ph", "iPr"):
        spec = CatalystSpec("Co", "bipyridine_pnnoh", substituent)
        for seed in (20270000, 20270817):
            for index, direction in enumerate(vectors):
                if time.time() >= args.deadline_epoch:
                    break
                name = f"{spec.catalyst_id}__seed{seed}__placement{index}"
                generated = generate_catalyst(spec, "active", seed=seed)
                metadata = generated.metadata()
                atoms = generated.atoms
                direction = np.asarray(direction, dtype=float)
                direction /= np.linalg.norm(direction)
                carbon, oxygen = metadata["carbonyl_indices"][0]
                atoms.positions[carbon] = 1.82 * direction
                atoms.positions[oxygen] = 2.96 * direction
                pairwise = atoms.get_all_distances() + np.eye(len(atoms)) * 1000
                if pairwise.min() < .75:
                    record = {"catalyst_id": spec.catalyst_id, "placement": index,
                              "seed": seed, "status": "starting_overlap_rejected", "utc": utc_now()}
                else:
                    folder = scratch / name
                    native = optimize_geometry(atoms, folder, charge=0, unpaired=0, threads=1,
                        solvent="toluene", convergence="vtight", max_cycles=500,
                        timeout_seconds=min(600, max(1, args.deadline_epoch - time.time())))
                    record = {"catalyst_id": spec.catalyst_id, "state": "active", "seed": seed,
                        "placement": index, "initial_carbonyl_unit_direction": direction.tolist(),
                        "conformer": name, "native_result": asdict(native), "metadata": metadata,
                        "solvent": "toluene", "solvation_state": "gsolv", "accepted_geometry_screen": False}
                    if native.converged:
                        optimized = read(folder / "xtbopt.xyz")
                        distances = optimized.get_distances(0, metadata["donor_indices"])
                        identity = catalyst_identity_screen(optimized, metadata)
                        retained = bool(np.all((distances > 1.4) & (distances < 3.1)))
                        record.update(chemical_identity=identity, coordination_retained=retained,
                            metal_donor_distances_A=distances.tolist(),
                            accepted_geometry_screen=retained and identity["retained"],
                            optimized_geometry_sha256=hashlib.sha256((folder / "xtbopt.xyz").read_bytes()).hexdigest())
                        if record["accepted_geometry_screen"]:
                            volume = buried_volume(optimized, 0, ligand_indices=metadata["ligand_indices"], n_samples=50000, seed=2026)
                            record.update(buried_volume_percent=volume.percent_buried_volume,
                                buried_volume_Wilson95_percent=list(volume.confidence_interval_percent),
                                terminal_donor_angle_deg=bite_angle(optimized, 0, metadata["donor_indices"][0], metadata["donor_indices"][2]),
                                donor_labels="P-N-N", P_M_P_angle_deg=None)
                            key = spec.catalyst_id + "|active"
                            if key not in best or native.energy_eV < best[key]["native_result"]["energy_eV"]:
                                best[key] = record
                    atomic_json(folder / "record.json", record)
                attempts += 1
                with (data / "records.jsonl").open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(record, allow_nan=False) + "\n")
                atomic_json(data / "best_conformers.json", best)
                print(json.dumps({"utc": utc_now(), "name": name,
                    "accepted": record.get("accepted_geometry_screen", False),
                    "identity": record.get("chemical_identity"), "error": record.get("status")}), flush=True)
    for key, record in best.items():
        name = key.replace("|", "__")
        atoms = read(Path(record["native_result"]["output_directory"]) / "xtbopt.xyz")
        result = characterize_stationary_point(atoms, data / "thermochemistry" / name,
            scratch / "thermochemistry" / name, charge=0, threads=1,
            deadline_epoch=args.deadline_epoch, solvent="toluene",
            coordination_indices=record["metadata"], identity_metadata=record["metadata"], source_metadata=record)
        print(json.dumps({"utc": utc_now(), "stationarity": name, "status": result["status"]}), flush=True)
    atomic_json(data / "summary.json", {"finished_utc": utc_now(), "attempts": attempts,
        "accepted_species": len(best), "method": "GFN2-xTB/ALPB toluene gsolv",
        "scope": "Alternative initial CO placements only; unrestricted optimized-identity gate unchanged"})


if __name__ == "__main__":
    main()
