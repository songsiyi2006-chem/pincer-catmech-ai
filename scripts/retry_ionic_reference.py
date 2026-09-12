"""Bounded real SCC warm-start recovery for a potassium tert-butoxide pair.

High electronic temperatures are diagnostic preparation only. Reported energy
and frequency calculations must return to 300 K electronic smearing.
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
from ase import Atoms
from ase.io import read
import numpy as np
from pincer_catmech.quantum.xtb_backend import _execute, executable_path, optimize_geometry, atomic_json, utc_now
from pincer_catmech.quantum.thermochemistry import characterize_stationary_point
from run_reference_thermochemistry import connectivity_check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deadline-epoch", required=True, type=float)
    args = parser.parse_args()
    source = REPO / "data/reference_thermochemistry/tert_butoxide_anion/stationarity/accepted_geometry.xyz"
    original = read(source)
    from rdkit import Chem
    graph = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)[O-]"))
    expected_bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in graph.GetBonds()]
    oxygen = next(i for i, atom in enumerate(original) if atom.symbol == "O")
    carbons = [i for i, atom in enumerate(original) if atom.symbol == "C"]
    carbon = min(carbons, key=lambda i: original.get_distance(i, oxygen))
    direction = original.positions[oxygen] - original.positions[carbon]
    direction /= np.linalg.norm(direction)
    out = REPO / "data/ionic_reference_retry"
    scratch = REPO.parent / "ionic-reference-scratch"
    out.mkdir(parents=True, exist_ok=True)
    result = {"started_utc": utc_now(), "source_geometry_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "method": "GFN2-xTB/ALPB toluene gsolv", "attempts": [], "accepted": False,
              "interpretation": "Single contact-pair model; not dissolved salt speciation or solubility"}
    for distance in (2.2, 2.6, 3.0, 3.5):
        atoms = original.copy()
        atoms += Atoms("K", positions=[atoms.positions[oxygen] + distance * direction])
        for temperature in (1000., 3000., 5000.):
            if time.time() >= args.deadline_epoch:
                break
            label = f"distance{distance:.1f}_warm{temperature:.0f}"
            folder = scratch / label
            item = {"initial_K_O_distance_A": distance, "warm_electronic_temperature_K": temperature,
                    "status": "started", "utc": utc_now()}
            try:
                code, output, wall, command = _execute(atoms, folder / "warm", executable_path(),
                    0, 0, 1, ["--grad"], min(120, args.deadline_epoch-time.time()), temperature, "toluene")
                item.update(warm_returncode=code, warm_wall_seconds=wall,
                            warm_output_sha256=hashlib.sha256(output.encode()).hexdigest())
                restart = folder / "warm/xtbrestart"
                if code or "normal termination of xtb" not in output.lower() or not restart.exists():
                    item["status"] = "warm_SCC_failed"
                else:
                    code, output, wall, command = _execute(atoms, folder / "cooled", executable_path(),
                        0, 0, 1, ["--grad"], min(120, args.deadline_epoch-time.time()), 300., "toluene", restart_file=restart)
                    item.update(cooled_returncode=code, cooled_wall_seconds=wall,
                                cooled_output_sha256=hashlib.sha256(output.encode()).hexdigest())
                    if code or "normal termination of xtb" not in output.lower():
                        item["status"] = "cooled_300K_SCC_failed"
                    else:
                        native = optimize_geometry(atoms, folder / "optimization", charge=0, unpaired=0,
                            threads=1, solvent="toluene", convergence="extreme", max_cycles=400,
                            timeout_seconds=min(600, args.deadline_epoch-time.time()),
                            restart_file=folder / "cooled/xtbrestart")
                        item["optimization"] = asdict(native)
                        item["status"] = native.status
                        if native.converged:
                            optimized = read(folder / "optimization/xtbopt.xyz")
                            stationarity = characterize_stationary_point(optimized, out / label,
                                folder / "stationarity", charge=0, threads=1, deadline_epoch=args.deadline_epoch,
                                solvent="toluene", source_metadata={"contact_pair": "tBuOK", "warm_start": item})
                            item["stationarity_status"] = stationarity["status"]
                            item["accepted"] = False
                            if stationarity["accepted"]:
                                final_atoms = read(out / label / "accepted_geometry.xyz")
                                item["chemical_identity"] = connectivity_check(final_atoms, expected_bonds)
                                item["final_K_O_distance_A"] = float(final_atoms.get_distance(oxygen, len(final_atoms)-1))
                                item["accepted"] = bool(item["chemical_identity"]["retained"] and item["final_K_O_distance_A"] < 4.0)
                            if item["accepted"]:
                                result.update(accepted=True, accepted_stationarity_file=str(out / label / "result.json"))
            except Exception as error:
                item.update(status="failed", error=f"{type(error).__name__}: {error}")
            result["attempts"].append(item)
            atomic_json(out / "summary.json", result)
            print(json.dumps(item, allow_nan=False), flush=True)
            if result["accepted"]:
                break
        if result["accepted"]:
            break
    result["finished_utc"] = utc_now()
    atomic_json(out / "summary.json", result)


if __name__ == "__main__":
    main()
