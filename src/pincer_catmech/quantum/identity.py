"""Geometric chemical-identity audit separate from mathematical stationarity.

No bond energies or oxidation-state assignments are inferred. Covalent contacts
are a documented radius screen, including every original ligand hydrogen.
"""
from functools import lru_cache
from typing import Mapping

from ase.data import covalent_radii
import numpy as np

from pincer_catmech.generators.combinatorial_pincer import CatalystSpec, _neutral_ligand


@lru_cache(maxsize=72)
def expected_ligand_graph(metal, backbone, substituent, state):
    """Recover generator atom order from the molecular graph without embedding."""
    from rdkit import Chem
    if state not in ("active", "hydrogenated", "protonated_reference"):
        raise ValueError("unknown catalyst state")
    molecule, maps, proton_map = _neutral_ligand(CatalystSpec(metal, backbone, substituent))
    if state == "active":
        site = maps[proton_map]
        proton = next(a.GetIdx() for a in molecule.GetAtomWithIdx(site).GetNeighbors() if a.GetAtomicNum() == 1)
        edited = Chem.RWMol(molecule)
        edited.RemoveAtom(proton)
        atom = edited.GetAtomWithIdx(site)
        atom.SetFormalCharge(-1)
        atom.SetNoImplicit(True)
        molecule = edited.GetMol()
        Chem.SanitizeMol(molecule)
    bonds = tuple((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1) for b in molecule.GetBonds())
    symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
    return bonds, symbols


def catalyst_identity_screen(atoms, metadata: Mapping, *, check_ancillaries=True):
    """Audit original molecular identity without constraining any calculation.

    Additional mapped substrate atoms in a reaction endpoint are excluded from
    the original-catalyst graph, allowing the specifically intended H uptake.
    Standalone hydrogenated-state metadata already includes its extra hydride
    and ligand proton, so their retention is checked explicitly.
    """
    bonds, symbols = expected_ligand_graph(*(metadata[key] for key in
        ("metal", "backbone", "substituent", "state")))
    ligand = list(metadata["ligand_indices"])
    if tuple(atoms[i].symbol for i in ligand) != symbols:
        raise ValueError("ligand atom order/elements differ from the generator graph")
    if metadata.get("ligand_bond_indices") is not None:
        supplied = {tuple(sorted(pair)) for pair in metadata["ligand_bond_indices"]}
        if supplied != {tuple(sorted(pair)) for pair in bonds}:
            raise ValueError("supplied ligand bonds conflict with the specified catalyst graph")
    expected = {tuple(sorted(pair)) for pair in bonds}
    indices = list(ligand)
    if check_ancillaries:
        for pair in metadata.get("carbonyl_indices", []):
            expected.add(tuple(sorted(pair)))
            indices.extend(pair)
        indices.extend(metadata.get("hydride_indices", []))
    if len(set(indices)) != len(indices):
        raise ValueError("original catalyst atom groups overlap")
    observed = set()
    for offset, i in enumerate(indices):
        for j in indices[:offset]:
            radius = covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number]
            if atoms.get_distance(i, j) < 1.28 * radius:
                observed.add(tuple(sorted((i, j))))
    lost, new = sorted(expected-observed), sorted(observed-expected)
    ligand_set = set(ligand)
    ligand_lost = [pair for pair in lost if set(pair) <= ligand_set]
    ligand_new = [pair for pair in new if set(pair) <= ligand_set]
    metal = metadata.get("metal_index", 0)
    hydrides = [float(atoms.get_distance(metal, i)) for i in metadata.get("hydride_indices", [])]
    carbonyls = [float(atoms.get_distance(metal, pair[0])) for pair in metadata.get("carbonyl_indices", [])]
    ancillary_retained = (all(1.0 < distance < 2.3 for distance in hydrides) and
                          all(1.3 < distance < 2.8 for distance in carbonyls)) if check_ancillaries else True
    return {"retained": not lost and not new and ancillary_retained,
            "ligand_retained": not ligand_lost and not ligand_new,
            "original_nonmetal_graph_retained": not lost and not new,
            "lost_bonds": lost, "new_bonds": new,
            "ligand_lost_bonds": ligand_lost, "ligand_new_bonds": ligand_new,
            "ancillary_coordination_retained": ancillary_retained,
            "metal_hydride_distances_A": hydrides, "metal_carbonyl_distances_A": carbonyls,
            "criterion": "Covalent graph at 1.28 times ASE covalent radii; original M-H 1.0-2.3 A and M-C(CO) 1.3-2.8 A",
            "interpretation": "Geometric intended-state screen; not bond-order, oxidation-state, spin or stability proof"}
