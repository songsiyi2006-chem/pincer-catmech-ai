"""Graph-mapped, directional explicit tBuOH docking around a hemiaminal.

The generated structures are unbiased optimization starting geometries, never
transition states. Open motifs are deliberate: a closed proton wire may require
a different substrate conformer and must be demonstrated by a TS/IRC calculation.
Coordinates are angstrom; energies in the association helpers are eV.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from ase import Atoms
from ase.data import covalent_radii, vdw_radii
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

HEMIAMINAL_SMILES = "OC(Nc1ccccc1)c1ccccc1"
TERT_BUTANOL_SMILES = "CC(C)(C)O"
EV_TO_KCAL_MOL = 23.060547830619


@dataclass
class SolvationCluster:
    atoms: Atoms
    atom_map: list[dict]
    fragments: list[list[int]]
    bonds: list[tuple[int, int]]
    hydrogen_bonds: list[dict]
    substrate_sites: dict[str, int]
    seed: int
    topology: str
    placement_trials: int

    def metadata(self) -> dict:
        closure_edges = []
        for fragment in self.fragments[1:]:
            for i, j in self.bonds:
                if i not in fragment or j not in fragment:
                    continue
                if {self.atoms[i].symbol, self.atoms[j].symbol} == {"O", "H"}:
                    donor, hydrogen = (i, j) if self.atoms[i].symbol == "O" else (j, i)
                    closure_edges += [{"donor": donor, "hydrogen": hydrogen, "acceptor": self.substrate_sites[key],
                                       "label": "unconstrained_closure_probe"} for key in ("oxygen", "nitrogen")]
        return {"atom_map": self.atom_map, "fragments": self.fragments,
                "bonds": self.bonds, "hydrogen_bonds": self.hydrogen_bonds,
                "substrate_sites": self.substrate_sites, "seed": self.seed,
                "topology": self.topology, "placement_trials": self.placement_trials,
                "charge": 0, "unpaired_electrons": 0,
                "coordinate_unit": "angstrom", "is_transition_state": False,
                "initial_hbond_metrics": hydrogen_bond_metrics(self.atoms, self.hydrogen_bonds),
                "initial_closure_probes": hydrogen_bond_metrics(self.atoms, closure_edges)}


def molecular_graph(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid molecular SMILES")
    return Chem.AddHs(mol)


def embed_molecule(smiles: str, seed: int = 20260913) -> Atoms:
    mol = molecular_graph(smiles)
    settings = AllChem.ETKDGv3()
    settings.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, settings) != 0:
        raise ValueError("RDKit embedding failed")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return Atoms([a.GetSymbol() for a in mol.GetAtoms()],
                 positions=mol.GetConformer().GetPositions())


def graph_bonds(mol: Chem.Mol) -> list[tuple[int, int]]:
    return [tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()]


def connectivity_check(atoms: Atoms, bonds) -> dict:
    expected = {tuple(sorted(pair)) for pair in bonds}
    observed = set()
    for i, j in combinations(range(len(atoms)), 2):
        if atoms.get_distance(i, j) < 1.28 * (covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]]):
            observed.add((i, j))
    return {"retained": observed == expected, "lost_bonds": sorted(expected - observed),
            "new_bonds": sorted(observed - expected),
            "criterion": "1.28 times ASE covalent-radius sum; all mapped atoms"}


def _validate_graph_geometry(atoms: Atoms, mol: Chem.Mol):
    if atoms.get_chemical_symbols() != [a.GetSymbol() for a in mol.GetAtoms()]:
        raise ValueError("geometry atom order does not match source graph")
    if np.any(atoms.pbc) or not np.all(np.isfinite(atoms.positions)):
        raise ValueError("finite isolated coordinates are required")
    if not connectivity_check(atoms, graph_bonds(mol))["retained"]:
        raise ValueError("geometry connectivity differs from mapped source graph")


def _donor(mol, element):
    candidates = [(a.GetIdx(), h.GetIdx()) for a in mol.GetAtoms()
                  if a.GetAtomicNum() == element for h in a.GetNeighbors()
                  if h.GetAtomicNum() == 1]
    if len(candidates) != 1:
        raise ValueError("expected exactly one mapped donor for requested element")
    return candidates[0]


def hydrogen_bond_metrics(atoms: Atoms, edges: list[dict]) -> list[dict]:
    records = []
    for edge in edges:
        donor, hydrogen, acceptor = (int(edge[k]) for k in ("donor", "hydrogen", "acceptor"))
        if atoms[hydrogen].symbol != "H" or atoms[acceptor].symbol not in ("O", "N"):
            raise ValueError("hydrogen-bond indices do not match chemical atom types")
        v = atoms.positions[donor] - atoms.positions[hydrogen]
        w = atoms.positions[acceptor] - atoms.positions[hydrogen]
        norm = np.linalg.norm(v) * np.linalg.norm(w)
        angle = float(np.degrees(np.arccos(np.clip(np.dot(v, w) / norm, -1, 1)))) if norm else 0.
        distance = float(np.linalg.norm(w))
        records.append({**edge, "H_acceptor_distance_A": distance,
                        "donor_H_acceptor_angle_deg": angle,
                        "within_initial_window": 1.6 - 1e-9 <= distance <= 2.1 + 1e-9 and angle >= 140.,
                        "retained_after_optimization": 1.4 <= distance <= 2.6 and angle >= 120.})
    return records


def steric_clearance(atoms: Atoms, fragment_ids, hydrogen_bonds) -> dict:
    """Reject interfragment overlaps; planned H...acceptor contacts are exempt.

    Other interfragment separations must exceed 0.65 times van der Waals radii
    sums. Full graph checks separately verify intrafragment covalent identity.
    """
    exempt = {tuple(sorted((e["hydrogen"], e["acceptor"]))) for e in hydrogen_bonds}
    i, j = np.triu_indices(len(atoms), 1)
    ids = np.asarray(fragment_ids)
    mask = ids[i] != ids[j]
    for a, b in exempt:
        mask &= ~((i == a) & (j == b))
    i, j = i[mask], j[mask]
    if not len(i):
        return {"accepted": True, "clashing_pairs": [], "closest_scaled_contact": None}
    threshold = .65 * (vdw_radii[atoms.numbers[i]] + vdw_radii[atoms.numbers[j]])
    distances = np.linalg.norm(atoms.positions[i] - atoms.positions[j], axis=1)
    ratios = distances / threshold
    index = int(np.argmin(ratios))
    worst = {"pair": [int(i[index]), int(j[index])], "distance_A": float(distances[index]),
             "threshold_A": float(threshold[index]), "ratio": float(ratios[index])}
    clashes = np.column_stack([i[ratios < 1], j[ratios < 1]]).tolist()
    return {"accepted": not clashes, "clashing_pairs": clashes, "closest_scaled_contact": worst}


def build_solvation_cluster(n_solvent: int, *, seed: int = 0,
                            substrate: Atoms | None = None,
                            solvent: Atoms | None = None,
                            hbond_distance_A: float = 1.85,
                            orientation_trials: int = 240) -> SolvationCluster:
    """Dock 1..3 explicit alcohols with NH/OH anchoring and a solvent relay.

    Seed parity changes the first anchor; with two alcohols both substrate
    donors are solvated; the third accepts the first alcohol's OH. Rotations
    search steric clearance, while a <=15 degree cone samples donor direction.
    This initializer uses no energy bias or assumed cyclic TS geometry.
    """
    if isinstance(n_solvent, bool) or not isinstance(n_solvent, int) or not 1 <= n_solvent <= 3:
        raise ValueError("n_solvent must be 1, 2 or 3")
    if not 1.6 <= hbond_distance_A <= 2.1 or orientation_trials < 1:
        raise ValueError("invalid hydrogen-bond distance or orientation trials")
    mol, alcohol = molecular_graph(HEMIAMINAL_SMILES), molecular_graph(TERT_BUTANOL_SMILES)
    substrate = substrate.copy() if substrate is not None else embed_molecule(HEMIAMINAL_SMILES)
    solvent = solvent.copy() if solvent is not None else embed_molecule(TERT_BUTANOL_SMILES)
    _validate_graph_geometry(substrate, mol)
    _validate_graph_geometry(solvent, alcohol)
    oh, nh = _donor(mol, 8), _donor(mol, 7)
    alcohol_o, alcohol_h = _donor(alcohol, 8)
    sites = {"oxygen": oh[0], "hydroxyl_hydrogen": oh[1], "nitrogen": nh[0], "amine_hydrogen": nh[1]}
    rng = np.random.default_rng(seed)
    total_trials = 0
    # Retry the whole placement when a locally good first alcohol blocks a relay.
    for attempt in range(12):
        atoms = substrate.copy()
        fragments = [list(range(len(atoms)))]
        fragment_ids = [0] * len(atoms)
        mapping = [{"global_index": i, "fragment": 0, "source_atom_index": i,
                    "symbol": a.symbol, "source_smiles": HEMIAMINAL_SMILES} for i, a in enumerate(atoms)]
        bonds = graph_bonds(mol)
        edges = []
        donors = [oh, nh] if seed % 2 == 0 else [nh, oh]
        success = True
        for number in range(n_solvent):
            donor, hydrogen = donors[number]
            unit = atoms.positions[hydrogen] - atoms.positions[donor]
            unit /= np.linalg.norm(unit)
            offset = len(atoms)
            edge = {"donor": donor, "hydrogen": hydrogen, "acceptor": offset + alcohol_o,
                    "label": "substrate_anchor" if number < 2 else "alcohol_relay"}
            best = None
            for trial in range(orientation_trials):
                total_trials += 1
                transverse = rng.normal(size=3)
                transverse -= unit * np.dot(transverse, unit)
                transverse /= max(np.linalg.norm(transverse), 1e-14)
                angle = np.deg2rad(rng.uniform(0, 15))
                direction = np.cos(angle) * unit + np.sin(angle) * transverse
                oxygen_position = atoms.positions[hydrogen] + hbond_distance_A * direction
                rotated = (solvent.positions - solvent.positions[alcohol_o]) @ Rotation.random(random_state=rng).as_matrix().T
                candidate = atoms + Atoms(solvent.symbols, positions=rotated + oxygen_position)
                clear = steric_clearance(candidate, fragment_ids + [number + 1] * len(solvent), edges + [edge])
                if not clear["accepted"]:
                    continue
                score = clear["closest_scaled_contact"]["ratio"]
                # Keep the first relay donor exposed for the three-solvent case.
                if n_solvent == 3 and number == 0:
                    probe = candidate.positions[offset + alcohol_h] + 1.85 * (
                        candidate.positions[offset + alcohol_h] - candidate.positions[offset + alcohol_o]) / np.linalg.norm(
                        candidate.positions[offset + alcohol_h] - candidate.positions[offset + alcohol_o])
                    score += .05 * min(np.linalg.norm(atoms.positions - probe, axis=1))
                if best is None or score > best[0]:
                    best = (score, candidate)
            if best is None:
                success = False
                break
            atoms = best[1]
            edges.append(edge)
            fragments.append(list(range(offset, len(atoms))))
            fragment_ids += [number + 1] * len(solvent)
            bonds += [(i + offset, j + offset) for i, j in graph_bonds(alcohol)]
            mapping += [{"global_index": offset + i, "fragment": number + 1, "source_atom_index": i,
                         "symbol": a.symbol, "source_smiles": TERT_BUTANOL_SMILES} for i, a in enumerate(solvent)]
            if number == 0:
                donors.append((offset + alcohol_o, offset + alcohol_h))
        if success and connectivity_check(atoms, bonds)["retained"]:
            return SolvationCluster(atoms, mapping, fragments, bonds, edges, sites, seed,
                "open_single_anchor" if n_solvent == 1 else "open_dual_anchor" if n_solvent == 2 else "open_dual_anchor_plus_alcohol_relay",
                total_trials)
    raise ValueError("no collision-free directional placement found within bounded search")


def association_energy(cluster_energy_eV, substrate_energy_eV, solvent_energies_eV) -> float:
    values = np.asarray([cluster_energy_eV, substrate_energy_eV, *solvent_energies_eV], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("association energies must be finite")
    return float(values[0] - np.sum(values[1:]))


def many_body_nonadditivity(cluster_energy_eV, monomer_energies_eV, pair_energies_eV) -> dict:
    """Frozen-geometry remainder beyond pair interactions (same protocol).

    ALPB fragment cavities are recalculated per subset; this includes nonadditive
    continuum response and must not be labeled pure gas-phase cooperativity.
    """
    monomers = np.asarray(monomer_energies_eV, dtype=float)
    expected = set(combinations(range(len(monomers)), 2))
    if set(pair_energies_eV) != expected:
        raise ValueError("all unique frozen fragment pairs are required")
    if not np.all(np.isfinite([cluster_energy_eV, *monomers, *pair_energies_eV.values()])):
        raise ValueError("finite frozen-fragment energies are required")
    interaction = float(cluster_energy_eV - monomers.sum())
    pairs = {pair: float(energy - monomers[pair[0]] - monomers[pair[1]])
             for pair, energy in pair_energies_eV.items()}
    return {"frozen_interaction_eV": interaction, "pair_sum_eV": sum(pairs.values()),
            "beyond_pair_nonadditivity_eV": interaction - sum(pairs.values()),
            "pair_interactions_eV": {f"{i}-{j}": e for (i, j), e in pairs.items()}}
