"""Actual ALPB/toluene GFN2-xTB reference minima for borrowing hydrogen.

RDKit only supplies starting coordinates. Every reported energy comes from
native xTB. Free ions and a potassium contact pair are separate reference
models; neither establishes dissolved t-BuOK speciation in toluene.
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

from ase import Atoms
from ase.data import covalent_radii
from ase.io import read, write
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pincer_catmech.quantum.xtb_backend import optimize_geometry, atomic_json, utc_now
from pincer_catmech.quantum.thermochemistry import characterize_stationary_point

SPECIES = {
    "benzyl_alcohol": "OCc1ccccc1",
    "benzaldehyde": "O=Cc1ccccc1",
    "aniline": "Nc1ccccc1",
    "hemiaminal": "OC(Nc1ccccc1)c1ccccc1",
    "imine": "c1ccc(/C=N/c2ccccc2)cc1",
    "product_amine": "c1ccc(CNc2ccccc2)cc1",
    "water": "O",
    "tert_butanol": "CC(C)(C)O",
    "tert_butoxide_anion": "CC(C)(C)[O-]",
    "potassium_tert_butoxide_contact_pair": "CC(C)(C)[O-]",
}


def assemble(smiles, seed, potassium_pair=False):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    settings = AllChem.ETKDGv3()
    settings.randomSeed = seed
    if AllChem.EmbedMolecule(molecule, settings) != 0:
        raise ValueError("RDKit embedding failed")
    if AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule, maxIters=500)
    atoms = Atoms([a.GetSymbol() for a in molecule.GetAtoms()],
                  positions=molecule.GetConformer().GetPositions())
    charge = Chem.GetFormalCharge(molecule)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()]
    if potassium_pair:
        oxygen = next(a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() == 8)
        carbon = molecule.GetAtomWithIdx(oxygen).GetNeighbors()[0].GetIdx()
        vector = atoms.positions[oxygen] - atoms.positions[carbon]
        vector /= np.linalg.norm(vector)
        atoms += Atoms("K", positions=[atoms.positions[oxygen] + 2.6 * vector])
        charge = 0
    return atoms, charge, bonds


def connectivity_check(atoms, bonds):
    # Include every covalent pair and forbid new covalent contacts. Potassium
    # coordination is recorded separately and is not a covalent topology test.
    expected = {tuple(sorted(pair)) for pair in bonds}
    observed = set()
    for i in range(len(atoms)):
        for j in range(i):
            if "K" in (atoms[i].symbol, atoms[j].symbol):
                continue
            radius_sum = covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number]
            if atoms.get_distance(i, j) < 1.28 * radius_sum:
                observed.add((j, i))
    return {"retained": expected == observed,
            "lost_bonds": sorted(expected - observed),
            "new_bonds": sorted(observed - expected),
            "criterion": "1.28 times ASE covalent-radius sum, excluding K coordination"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deadline-epoch", type=float, required=True)
    parser.add_argument("--conformers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if not 1 <= args.conformers <= 50 or not 1 <= args.threads <= 8:
        raise ValueError("invalid resource limits")
    out = REPO / "data" / "reference_thermochemistry"
    scratch = REPO.parent / "reference-scratch"
    out.mkdir(parents=True, exist_ok=True)
    summary = {"started_utc": utc_now(), "method": "GFN2-xTB", "solvent": "toluene",
               "solvation_state": "gsolv", "species": {},
               "scope": "Molecular reference minima; ion-pair solubility and speciation are not determined"}
    for name, smiles in SPECIES.items():
        if time.time() >= args.deadline_epoch:
            break
        species_dir = out / name
        species_dir.mkdir(exist_ok=True)
        if (species_dir / "summary.json").exists():
            summary["species"][name] = json.loads((species_dir / "summary.json").read_text())
            continue
        records, best = [], None
        for index in range(args.conformers):
            if time.time() >= args.deadline_epoch:
                break
            atoms, charge, bonds = assemble(smiles, 20260 + index, "contact_pair" in name)
            folder = scratch / name / f"conformer_{index:02d}"
            native = optimize_geometry(atoms, folder, charge=charge, unpaired=0,
                threads=args.threads, max_cycles=400, convergence="extreme", solvent="toluene",
                timeout_seconds=min(600, max(1, args.deadline_epoch - time.time())))
            record = {"conformer": index, "smiles": smiles, "native_result": asdict(native),
                      "charge": charge, "unpaired_electrons": 0}
            if native.converged:
                optimized = read(folder / "xtbopt.xyz")
                record["connectivity"] = connectivity_check(optimized, bonds)
                if record["connectivity"]["retained"] and (best is None or native.energy_eV < best[0]):
                    best = (native.energy_eV, optimized, charge, bonds, folder, index)
            records.append(record)
            print(json.dumps({"utc": utc_now(), "species": name, "optimization": record}), flush=True)
        item = {"species": name, "smiles": smiles, "conformers": records,
                "status": "no_converged_reference", "accepted": False}
        if best is not None:
            _, atoms, charge, bonds, folder, index = best
            write(species_dir / "best_optimized.xyz", atoms)
            item.update(best_conformer=index, best_geometry_sha256=hashlib.sha256(
                (species_dir / "best_optimized.xyz").read_bytes()).hexdigest())
            native = characterize_stationary_point(atoms, species_dir / "stationarity",
                scratch / name / "stationarity", charge=charge, threads=args.threads,
                max_repairs=2, deadline_epoch=args.deadline_epoch, solvent="toluene",
                source_metadata={"species": name, "smiles": smiles,
                    "optimization_output": str(folder / "xtb.out"),
                    "state_model": "contact ion pair" if "contact_pair" in name else "isolated molecular species"})
            item.update(status=native["status"], accepted=native["accepted"],
                        thermochemistry_file=str(species_dir / "stationarity" / "result.json"))
            if native["accepted"]:
                retained = connectivity_check(read(species_dir / "stationarity" / "accepted_geometry.xyz"), bonds)
                item["final_connectivity"] = retained
                if not retained["retained"]:
                    item.update(status="reference_identity_changed", accepted=False)
        atomic_json(species_dir / "summary.json", item)
        summary["species"][name] = item
        atomic_json(out / "summary.json", summary)
        print(json.dumps({"utc": utc_now(), "species": name, "status": item["status"]}), flush=True)
    summary.update(finished_utc=utc_now(), accepted_species=sum(x["accepted"] for x in summary["species"].values()))
    atomic_json(out / "summary.json", summary)


if __name__ == "__main__":
    main()
