"""One-worker actual GFN2-xTB ALPB condensation search, with atomic provenance.

The first model is a neutral hemiaminal/tert-butanol proton shuttle, not a
base-concentration model. In the optional anionic model tBuO- accepts the N-H
proton while a separately mapped tBuOH donates its proton to the leaving O.
The two shuttle carbon skeletons exchange protonation state; charge stays -1.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/"src"))

import numpy as np
from ase.io import read, write
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from rdkit import Chem

from pincer_catmech.kinetics.condensation_search import run_condensation_path, observed_bonds
from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now


def reference(name, smiles):
    folder = REPO/"data"/"reference_thermochemistry"/name/"stationarity"
    outcome = json.loads((folder/"result.json").read_text(encoding="utf-8"))
    if outcome.get("status") != "verified_minimum" or outcome.get("accepted") is not True or outcome.get("solvent") != "toluene":
        raise ValueError("Reference lacks accepted ALPB minimum evidence: "+name)
    path = folder/"accepted_geometry.xyz"
    atoms = read(path)
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if [a.GetAtomicNum() for a in molecule.GetAtoms()] != atoms.numbers.tolist():
        raise ValueError("Reference atom ordering differs from supplied molecular graph")
    bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in molecule.GetBonds()}
    if observed_bonds(atoms) != bonds:
        raise ValueError("Reference geometry and annotated covalent graph differ")
    return atoms, molecule, bonds, {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def assemble(seed=2026, anionic=False):
    hemi, hm, hb, hsource = reference("hemiaminal", "OC(Nc1ccccc1)c1ccccc1")
    alcohol, am, ab, asource = reference("tert_butanol", "CC(C)(C)O")
    oxygen, carbon, nitrogen = 0, 1, 2
    hn = next(a.GetIdx() for a in hm.GetAtomWithIdx(nitrogen).GetNeighbors() if a.GetAtomicNum() == 1)
    ho = next(a.GetIdx() for a in hm.GetAtomWithIdx(oxygen).GetNeighbors() if a.GetAtomicNum() == 1)
    ao = next(a.GetIdx() for a in am.GetAtoms() if a.GetAtomicNum() == 8)
    ah = next(a.GetIdx() for a in am.GetAtomWithIdx(ao).GetNeighbors() if a.GetAtomicNum() == 1)
    local = alcohol.positions-alcohol.positions[ao]
    rng = np.random.default_rng(seed)
    targets = hemi.positions[[hn, oxygen]]

    def positions(parameters):
        return local@Rotation.from_rotvec(parameters[:3]).as_matrix().T+parameters[3:]

    def objective(parameters):
        placed = positions(parameters)
        distances = np.linalg.norm(hemi.positions[:, None]-placed[None], axis=2)
        # Geometric seed score only; no energy is inferred from this function.
        clash = np.maximum(1.45-distances, 0)
        return float(150*np.sum(clash**2)+(0 if anionic else 30)*(np.linalg.norm(placed[ao]-targets[0])-1.75)**2
                     +30*(np.linalg.norm(placed[ah]-targets[1])-1.75)**2
                     +(0 if anionic else 5)*(np.linalg.norm(placed[ao]-hemi.positions[nitrogen])-2.8)**2)
    choices = []
    center = targets.mean(axis=0)
    for _ in range(24):
        initial = np.r_[Rotation.random(random_state=rng).as_rotvec(), center+rng.normal(size=3)]
        fitted = minimize(objective, initial, method="BFGS", options={"maxiter": 180})
        choices.append((float(fitted.fun), fitted.x))
    score, best = min(choices, key=lambda row: row[0])
    alcohol.positions = positions(best)
    offset = len(hemi)
    ot, ht = offset+ao, offset+ah
    reactant = hemi+alcohol
    rbonds = hb|{(a+offset, b+offset) for a, b in ab}
    metadata = {"created_utc": utc_now(), "seed": seed, "sources": [hsource, asource],
        "assembly_geometric_score_only": score, "potassium_explicit": False,
        "interpretation": "Neutral tBuOH proton shuttle; no base kinetic order or tBuOK speciation inferred"}
    receiver = ot
    if anionic:
        base, bm, bb, bsource = reference("tert_butoxide_anion", "CC(C)(C)[O-]")
        bo = next(a.GetIdx() for a in bm.GetAtoms() if a.GetAtomicNum() == 8)
        base_local = base.positions-base.positions[bo]
        extension = hemi.positions[hn]-hemi.positions[nitrogen]
        base_target = hemi.positions[hn]+1.70*extension/np.linalg.norm(extension)
        def base_positions(parameters):
            return base_local@Rotation.from_rotvec(parameters[:3]).as_matrix().T+parameters[3:]
        def base_objective(parameters):
            placed = base_positions(parameters)
            distances = np.linalg.norm(reactant.positions[:, None]-placed[None], axis=2)
            return float(150*np.sum(np.maximum(1.45-distances, 0)**2)
                         +50*np.sum((placed[bo]-base_target)**2))
        base_choices = []
        for _ in range(24):
            fitted = minimize(base_objective, np.r_[Rotation.random(random_state=rng).as_rotvec(), base_target],
                              method="BFGS", options={"maxiter": 180})
            base_choices.append((float(fitted.fun), fitted.x))
        base_score, base_parameters = min(base_choices, key=lambda row: row[0])
        base.positions = base_positions(base_parameters)
        boffset = len(reactant)
        receiver = boffset+bo
        reactant += base
        rbonds |= {(a+boffset, b+boffset) for a, b in bb}
        metadata["sources"].append(bsource)
        metadata.update(base_assembly_geometric_score_only=base_score,
            interpretation="Isolated anionic tBuO-/tBuOH cluster without potassium; N-H protonates the original tBuO- while original tBuOH protonates leaving O; mapped shuttle carbon skeletons exchange protonation states")
    product = reactant.copy()
    direction = reactant.positions[oxygen]-reactant.positions[carbon]
    direction /= np.linalg.norm(direction)
    shift = (3.1-reactant.get_distance(oxygen, carbon))*direction
    product.positions[[oxygen, ho]] += shift
    vector = reactant.positions[nitrogen]-reactant.positions[receiver]
    product.positions[hn] = product.positions[receiver]+1.00*vector/np.linalg.norm(vector)
    vector = product.positions[ot]-product.positions[oxygen]
    product.positions[ht] = product.positions[oxygen]+0.98*vector/np.linalg.norm(vector)
    pbonds = set(rbonds)
    for pair in ((oxygen, carbon), (nitrogen, hn), (ot, ht)):
        pbonds.remove(tuple(sorted(pair)))
    for pair in ((receiver, hn), (oxygen, ht)):
        pbonds.add(tuple(sorted(pair)))
    coordinates = {
        "N_H_to_acceptor_O": {"terms": [(nitrogen, hn, 1), (receiver, hn, -1)], "minimum_overlap": 0.10},
        "shuttle_H_to_leaving_O": {"terms": [(ot, ht, 1), (oxygen, ht, -1)], "minimum_overlap": 0.10},
        "C_O_cleavage": {"terms": [(carbon, oxygen, 1)], "minimum_overlap": 0.05},
        "C_N_contraction": {"terms": [(carbon, nitrogen, -1)], "minimum_overlap": 0.03}}
    metadata.update(mapped_indices={"C": carbon, "N": nitrogen, "leaving_O": oxygen,
        "N_H": hn, "original_O_H": ho, "proton_acceptor_O": receiver, "shuttle_O": ot, "shuttle_H": ht},
        total_atoms=len(reactant), total_charge=-1 if anionic else 0)
    return reactant, product, [sorted(rbonds), sorted(pbonds)], coordinates, (carbon, nitrogen), metadata


def product_from_reference_components(seed_atoms, *, anionic, water_distance):
    """Fit actual separately accepted product geometries to the mapped seed.

    Bond-order-correct imine coordinates avoid seeding a product with the
    reactant's C-N single-bond geometry. Rigid fitting is geometry preparation,
    never an energy evaluation or a substitute for complete-cluster optimization.
    """
    pieces = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for smiles in
              ("OC(Nc1ccccc1)c1ccccc1", "CC(C)(C)O")]
    if anionic:
        pieces.append(Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)[O-]")))
    graph = pieces[0]
    for piece in pieces[1:]:
        graph = Chem.CombineMols(graph, piece)
    if [atom.GetAtomicNum() for atom in graph.GetAtoms()] != seed_atoms.numbers.tolist():
        raise ValueError("Product reference mapping has incompatible atomic identities")
    rw = Chem.RWMol(graph)
    nitrogen, carbon, oxygen = 2, 1, 0
    hn = next(a.GetIdx() for a in rw.GetAtomWithIdx(nitrogen).GetNeighbors() if a.GetAtomicNum() == 1)
    ot = pieces[0].GetNumAtoms()+next(a.GetIdx() for a in pieces[1].GetAtoms() if a.GetAtomicNum() == 8)
    ht = next(a.GetIdx() for a in rw.GetAtomWithIdx(ot).GetNeighbors() if a.GetAtomicNum() == 1)
    receiver = ot if not anionic else pieces[0].GetNumAtoms()+pieces[1].GetNumAtoms()+next(a.GetIdx() for a in pieces[2].GetAtoms() if a.GetAtomicNum() == 8)
    for a, b in ((oxygen, carbon), (nitrogen, hn), (ot, ht)):
        rw.RemoveBond(a, b)
    rw.AddBond(receiver, hn, Chem.BondType.SINGLE)
    rw.AddBond(oxygen, ht, Chem.BondType.SINGLE)
    rw.GetBondBetweenAtoms(carbon, nitrogen).SetBondType(Chem.BondType.DOUBLE)
    if anionic:
        rw.GetAtomWithIdx(receiver).SetFormalCharge(0)
        rw.GetAtomWithIdx(ot).SetFormalCharge(-1)
    molecule = rw.GetMol()
    Chem.SanitizeMol(molecule)
    mappings = []
    fragments = Chem.GetMolFrags(molecule, asMols=True, fragsMolAtomMapping=mappings)
    references = [reference(name, smiles) for name, smiles in (
        ("imine", "c1ccc(/C=N/c2ccccc2)cc1"), ("water", "O"),
        ("tert_butanol", "CC(C)(C)O"), ("tert_butoxide_anion", "CC(C)(C)[O-]"))]
    result = seed_atoms.copy()
    used = []
    for fragment, indices in zip(fragments, mappings):
        target = seed_atoms.positions[list(indices)]
        fits = []
        for source_atoms, source_mol, _, provenance in references:
            if len(source_atoms) != len(indices):
                continue
            for match in source_mol.GetSubstructMatches(fragment, useChirality=False, uniquify=False, maxMatches=256):
                if len(match) != len(indices):
                    continue
                source = source_atoms.positions[list(match)]
                centered = source-source.mean(axis=0)
                u, _, vt = np.linalg.svd(centered.T@(target-target.mean(axis=0)))
                rotation = u@vt
                if np.linalg.det(rotation) < 0:
                    u[:, -1] *= -1
                    rotation = u@vt
                positioned = centered@rotation+target.mean(axis=0)
                fits.append((float(np.sum((positioned-target)**2)), positioned, provenance, match))
        if not fits:
            raise ValueError("No atom-mapped accepted reference matches a product component")
        _, positioned, provenance, match = min(fits, key=lambda item: item[0])
        result.positions[list(indices)] = positioned
        used.append({**provenance, "target_indices": list(indices), "source_indices": list(match)})
        if oxygen in indices:
            water_indices = list(indices)
    direction = result.positions[oxygen]-result.positions[carbon]
    result.positions[water_indices] += (water_distance-np.linalg.norm(direction))*direction/np.linalg.norm(direction)
    return result, used


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb", required=True)
    parser.add_argument("--deadline-epoch", type=float, required=True)
    parser.add_argument("--budget-seconds", type=float, default=1800)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--anionic", action="store_true")
    parser.add_argument("--restart-search", type=Path)
    parser.add_argument("--restart-step", type=int)
    parser.add_argument("--activated-from-search", type=Path)
    parser.add_argument("--certify-endpoints", action="store_true")
    parser.add_argument("--product-reference-geometry", action="store_true")
    parser.add_argument("--product-water-distance", type=float, default=4.3)
    parser.add_argument("--neb-steps", type=int, default=180)
    parser.add_argument("--neb-maxstep", type=float, default=0.10)
    args = parser.parse_args()
    if args.budget_seconds <= 0:
        raise ValueError("budget must be positive")
    if args.activated_from_search and (not args.anionic or args.restart_search):
        raise ValueError("Activated cluster requires --anionic and cannot also restart a band")
    tag = ("activated_anionic_elimination" if args.activated_from_search else
           "anionic_tBuO_tBuOH" if args.anionic else "neutral_tBuOH")+f"_seed{args.seed}"
    directory = REPO/"data"/"condensation_search"/tag
    scratch = REPO.parent/"condensation-scratch"/tag
    if directory.exists() or scratch.exists():
        raise FileExistsError("Use a new seed to preserve this search's immutable evidence")
    started = time.time()
    restart_band = None
    selected_restart_frames = None
    if args.activated_from_search is not None:
        prior = args.activated_from_search.resolve()
        protocol = json.loads((prior/"protocol.json").read_text(encoding="utf-8"))
        previous = json.loads((prior/"result.json").read_text(encoding="utf-8"))
        if protocol["charge"] != -1 or protocol["solvent"] != "toluene" or previous["status"] == "running":
            raise ValueError("Require completed actual anionic ALPB search evidence")
        reactant_path = prior/"endpoint_0.xyz"
        reactant = read(reactant_path)
        old_assembly = prior.parent/(prior.name+"_assembly")
        product_path = old_assembly/"product_seed.xyz"
        product = read(product_path)
        coordinates = dict(protocol["reaction_coordinates"])
        transfer = coordinates.pop("N_H_to_acceptor_O")
        donor, hydrogen, plus = transfer["terms"][0]
        acceptor, same_hydrogen, minus = transfer["terms"][1]
        if hydrogen != same_hydrogen or plus != 1 or minus != -1:
            raise ValueError("Previous annotated proton-transfer identity is not supported")
        activated_bonds = {tuple(sorted(pair)) for pair in protocol["endpoint_bonds"][0]}
        activated_bonds.remove(tuple(sorted((donor, hydrogen))))
        activated_bonds.add(tuple(sorted((acceptor, hydrogen))))
        if observed_bonds(reactant) != activated_bonds:
            raise ValueError("Actual endpoint changed beyond the explicitly specified N-H acid-base transfer")
        bonds = [sorted(activated_bonds), protocol["endpoint_bonds"][1]]
        cn = protocol["C_N_indices"]
        metadata = {"created_utc": utc_now(), "total_charge": -1, "potassium_explicit": False,
            "interpretation": "Isolated anionic elimination from N-deprotonated hemiaminal plus two mapped tBuOH molecules to imine/water/tBuOH/tBuO-. Requires a separate preceding acid-base equilibrium; not the original neutral-hemiaminal elementary step or a total-base rate law.",
            "actual_activated_source": str(reactant_path),
            "actual_activated_source_sha256": hashlib.sha256(reactant_path.read_bytes()).hexdigest(),
            "product_seed_source": str(product_path),
            "product_seed_source_sha256": hashlib.sha256(product_path.read_bytes()).hexdigest(),
            "explicit_proton_transfer": {"donor": donor, "hydrogen": hydrogen, "acceptor": acceptor}}
    elif args.restart_search is None:
        reactant, product, bonds, coordinates, cn, metadata = assemble(args.seed, args.anionic)
    else:
        prior = args.restart_search.resolve()
        protocol = json.loads((prior/"protocol.json").read_text(encoding="utf-8"))
        outcome = json.loads((prior/"result.json").read_text(encoding="utf-8"))
        if outcome["status"] == "running":
            raise ValueError("Do not start a second worker from a still-running search")
        if protocol["charge"] != (-1 if args.anionic else 0) or protocol["solvent"] != "toluene":
            raise ValueError("Restart charge or solvent convention differs")
        restart_band = prior/"final_band.xyz"
        if not restart_band.exists():
            restart_band = prior/"current_band.xyz"
        if args.restart_step is None:
            frames = read(restart_band, index=":")
        else:
            from ase.io.trajectory import Trajectory
            trajectory_path = prior/"neb.traj"
            with Trajectory(trajectory_path) as trajectory:
                first = 9*args.restart_step
                if args.restart_step < 0 or first+9 > len(trajectory):
                    raise ValueError("Restart step lies outside the saved complete nine-frame bands")
                frames = [trajectory[i].copy() for i in range(first, first+9)]
            selected_restart_frames = frames
        reactant, product = frames[0], frames[-1]
        bonds, coordinates, cn = protocol["endpoint_bonds"], protocol["reaction_coordinates"], protocol["C_N_indices"]
        metadata = {"created_utc": utc_now(), "interpretation": protocol["model"],
            "restart_protocol": str(prior/"protocol.json"),
            "restart_protocol_sha256": hashlib.sha256((prior/"protocol.json").read_bytes()).hexdigest(),
            "previous_outcome": outcome["status"], "new_seed_geometry_generated": False}
        if selected_restart_frames is not None:
            metadata.update(restart_trajectory=str(trajectory_path), restart_step=args.restart_step,
                restart_trajectory_sha256=hashlib.sha256(trajectory_path.read_bytes()).hexdigest(),
                restart_frame_indices=list(range(first, first+9)))
    preparation = directory.parent/(tag+"_assembly")
    preparation.mkdir(parents=True, exist_ok=False)
    if args.product_reference_geometry:
        if args.restart_search:
            raise ValueError("Reference replacement cannot alter endpoints of a restarted band")
        if not np.isfinite(args.product_water_distance) or args.product_water_distance <= 2.5:
            raise ValueError("Product water distance must be finite and above 2.5 angstrom")
        product, provenance = product_from_reference_components(product, anionic=args.anionic,
            water_distance=args.product_water_distance)
        metadata.update(product_component_references=provenance,
                        product_water_carbon_seed_distance_angstrom=args.product_water_distance)
    if selected_restart_frames is not None:
        restart_band = preparation/f"restart_step_{args.restart_step}.xyz"
        write(restart_band, selected_restart_frames)
    atomic_json(preparation/"metadata.json", metadata)
    write(preparation/"reactant_seed.xyz", reactant)
    write(preparation/"product_seed.xyz", product)
    result = run_condensation_path(reactant, product, output_dir=directory, scratch_dir=scratch,
        executable=args.xtb, charge=-1 if args.anionic else 0, endpoint_bonds=bonds,
        reaction_coordinates=coordinates, cn=cn, deadline_epoch=min(args.deadline_epoch, started+args.budget_seconds),
        model_label=metadata["interpretation"], restart_band=restart_band, neb_steps=args.neb_steps,
        neb_maxstep=args.neb_maxstep, certify_endpoints=args.certify_endpoints)
    print(json.dumps({"status": result["status"], "accepted": result["accepted"],
                      "failure_reasons": result["failure_reasons"], "directory": str(directory)}), flush=True)


if __name__ == "__main__":
    main()
