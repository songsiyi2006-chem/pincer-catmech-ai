"""Export reaction thermodynamics from accepted actual ALPB reference minima.

This postprocessing script does not run quantum chemistry, fit a model, or
produce activation barriers. All source hashes and spectra are revalidated.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
import csv
import hashlib
import json
import math
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/"src"))
import numpy as np
from ase.data import covalent_radii
from ase.io import read
from rdkit import Chem
from pincer_catmech.quantum.xtb_backend import atomic_json

TEMPERATURES = (298.15, 340., 350., 360., 370., 380., 383.15, 390., 400., 410., 420., 430., 440.)
REACTIONS = {
    "hemiaminal_formation": {"hemiaminal": 1, "benzaldehyde": -1, "aniline": -1},
    "hemiaminal_dehydration": {"imine": 1, "water": 1, "hemiaminal": -1},
    "net_borrowing_hydrogen": {"product_amine": 1, "water": 1, "benzyl_alcohol": -1, "aniline": -1},
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validated_source(name):
    folder = REPO/"data"/"reference_thermochemistry"/name/"stationarity"
    result_path = folder/"result.json"
    record = json.loads(result_path.read_text(encoding="utf-8"))
    if record.get("status") != "verified_minimum" or record.get("accepted") is not True:
        raise ValueError(name+" is not an accepted minimum")
    if record.get("method") != "GFN2-xTB" or record.get("solvent") != "toluene" or record.get("solvation_state") != "gsolv":
        raise ValueError(name+" has an inconsistent calculation protocol")
    geometry_path = folder/"accepted_geometry.xyz"
    if sha(geometry_path) != record["accepted_geometry_sha256"]:
        raise ValueError(name+" accepted geometry hash changed")
    atoms = read(geometry_path)
    molecule = Chem.AddHs(Chem.MolFromSmiles(record["source_metadata"]["smiles"]))
    if [atom.GetAtomicNum() for atom in molecule.GetAtoms()] != atoms.numbers.tolist():
        raise ValueError(name+" source graph and geometry have different atomic identities")
    expected_bonds = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in molecule.GetBonds()}
    observed = {(j, i) for i in range(len(atoms)) for j in range(i) if atoms.get_distance(i, j) <
                1.28*(covalent_radii[atoms.numbers[i]]+covalent_radii[atoms.numbers[j]])}
    if expected_bonds != observed:
        raise ValueError(name+" minimum does not retain its annotated molecular graph")
    table = {float(row["temperature"]): row for row in record["thermochemistry"]}
    if len(table) != len(record["thermochemistry"]) or tuple(sorted(table)) != TEMPERATURES:
        raise ValueError(name+" has an incomplete or duplicate temperature table")
    sources = {row["source"] for row in table.values()}
    if len(sources) != 1:
        raise ValueError(name+" temperature table mixes Hessian geometries")
    original_hessian_source = next(iter(sources))
    # Sources preserve the original workstation path. The certified NPZ is
    # shipped beside its result; the exact hash below also checks relocation.
    hessian_path = folder/PureWindowsPath(original_hessian_source).name
    hessian_hash = sha(hessian_path)
    with np.load(hessian_path) as hessian:
        frequencies = hessian["frequencies_cm1"]
        if not np.array_equal(atoms.numbers, hessian["atomic_numbers"]) or not np.allclose(atoms.positions, hessian["positions"], atol=1e-6, rtol=0):
            raise ValueError(name+" accepted geometry and Hessian geometry differ")
        if len(frequencies) != 3*len(atoms)-int(hessian["external_rank"]) or np.any(frequencies <= 0):
            raise ValueError(name+" lacks a complete all-positive internal spectrum")
        if float(hessian["antisymmetry_relative"]) > 0.05:
            raise ValueError(name+" Hessian exceeds the numerical consistency threshold")
    for temperature, row in table.items():
        if row["source_sha256"] != hessian_hash or not np.allclose(row["frequencies_cm1"], frequencies, rtol=0, atol=1e-8):
            raise ValueError(name+" reported source hash or spectrum differs from its NPZ")
        if row["concentration_mol_l"] != 1 or row["pressure_pa"] != 101325 or row["rotor_inertia"] != "grimme" or row["cutoff_cm1"] != 100:
            raise ValueError(name+" uses an inconsistent thermochemistry convention")
        reconstructed = row["H_298"]-temperature*row["S_qRRHO"]/1000+row["delta_G_concentration"]
        if not math.isfinite(row["G_298_qRRHO_sol"]) or not math.isclose(reconstructed, row["G_298_qRRHO_sol"], rel_tol=0, abs_tol=1e-8):
            raise ValueError(name+" Gibbs-energy bookkeeping does not reconstruct")
    entry = {"species_id": name, "formula": atoms.get_chemical_formula(), "charge": record["charge"],
        "smiles": record["source_metadata"]["smiles"],
        "result_path": str(result_path.relative_to(REPO)), "result_sha256": sha(result_path),
        "geometry_path": str(geometry_path.relative_to(REPO)), "geometry_sha256": sha(geometry_path),
        "hessian_path": str(hessian_path.relative_to(REPO)), "hessian_sha256": hessian_hash,
        "original_hessian_source": original_hessian_source,
        "internal_mode_count": len(frequencies), "minimum_frequency_cm1": float(frequencies.min()),
        "accepted_minimum": True, "molecular_identity_checked": True,
        "metadata_correction_audit": record.get("metadata_correction_audit"),
        "temperature_field": "temperature", "G_field": "G_298_qRRHO_sol",
        "G_field_interpretation": "Model G at the row temperature; legacy 298 name is not a temperature override"}
    return table, entry, Counter(atoms.numbers.tolist())


def main():
    names = sorted({name for reaction in REACTIONS.values() for name in reaction})
    validated = {name: validated_source(name) for name in names}
    for label, stoichiometry in REACTIONS.items():
        elements = {z for name in stoichiometry for z in validated[name][2]}
        if any(sum(coefficient*validated[name][2][z] for name, coefficient in stoichiometry.items()) for z in elements):
            raise ValueError(label+" violates atom balance")
        if sum(coefficient*validated[name][1]["charge"] for name, coefficient in stoichiometry.items()):
            raise ValueError(label+" violates charge balance")
    rows = [{"temperature_K": temperature, **{
        label+"_delta_G_kcal_mol": math.fsum(coefficient*validated[name][0][temperature]["G_298_qRRHO_sol"]
            for name, coefficient in stoichiometry.items()) for label, stoichiometry in REACTIONS.items()}}
        for temperature in TEMPERATURES]
    directory = REPO/"data"/"datasets"
    directory.mkdir(parents=True, exist_ok=True)
    target = directory/"reference_reaction_thermodynamics.csv"
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(), "status": "verified_postprocessing",
        "csv_path": str(target.relative_to(REPO)), "csv_sha256": sha(target),
        "generator": str(Path(__file__).relative_to(REPO)), "generator_sha256": sha(__file__),
        "temperatures_K": TEMPERATURES, "row_count": len(rows), "energy_unit": "kcal/mol",
        "method": "GFN2-xTB + ALPB toluene/gsolv + entropy-only Grimme qRRHO + 1 M standard state",
        "reaction_stoichiometry": REACTIONS,
        "reaction_delta_n": {label: sum(stoichiometry.values()) for label, stoichiometry in REACTIONS.items()},
        "definition": "delta_G_standard = sum(products G_standard) - sum(reactants G_standard)",
        "sources": [validated[name][1] for name in names],
        "interpretation": "Reaction thermodynamics of accepted independently optimized molecular references. No transition state, activation barrier, rate constant, TOF, or equilibrium population with ion-pair speciation is inferred.",
        "temperature_approximation": "Fixed accepted geometry/Hessian and fixed ALPB excess contribution reused; no independent solvation enthalpy or solvent temperature derivative",
        "standard_state": "ALPB gsolv excess contribution is already in the potential; one ideal 1 atm to 1 M shift per molecular species at each actual T",
        "checks": ["accepted minimum flags", "exact XYZ hash", "exact NPZ hash", "XYZ/NPZ positions and atom order",
                   "all-positive complete internal spectrum", "every row spectrum equals NPZ", "graph identity",
                   "consistent 13T and reference-state conventions", "G reconstruction", "atom and charge balances"]}
    atomic_json(directory/"reference_reaction_thermodynamics.manifest.json", manifest)
    print(json.dumps(next(row for row in rows if row["temperature_K"] == 383.15)))


if __name__ == "__main__":
    main()
