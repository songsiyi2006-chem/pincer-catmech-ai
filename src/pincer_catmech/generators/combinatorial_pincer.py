"""RDKit constrained organic-ligand assembly into explicit metal complexes.

These are designed starting structures, not experimentally authenticated catalysts
or optimized minima. No electronic energies are fabricated or returned. The
bipyridine-OH scaffold follows de Boer et al., Organometallics 2017,
doi:10.1021/acs.organomet.7b00111; its Ph/iPr and cross-metal analogues here are
proposals. Its responsive sites are OH and benzylic CH, not two NH groups.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Literal

import numpy as np
from ase import Atoms
from ase.io import write

BACKBONES = ("macho_pnp", "pyridine_pnn", "bipyridine_pnnoh")
METALS = ("Ru", "Mn", "Fe", "Co")
SUBSTITUENTS = ("Ph", "iPr")
STATES = ("active", "hydrogenated", "protonated_reference")
OXIDATION_STATES = {"Ru": 2, "Mn": 1, "Fe": 2, "Co": 1}
BACKBONE_REFERENCES = {
    "macho_pnp": "MACHO-type HN(CH2CH2PR2)2 connectivity; designed cross-metal ancillary set",
    "pyridine_pnn": "Pyridine PNN with phosphinomethyl and diethylaminomethyl arms; deprotonated benzylic-carbon resonance model",
    "bipyridine_pnnoh": "https://doi.org/10.1021/acs.organomet.7b00111",
}


@dataclass(frozen=True)
class CatalystSpec:
    metal: str
    backbone: str
    substituent: str

    def __post_init__(self):
        if self.metal not in METALS or self.backbone not in BACKBONES or self.substituent not in SUBSTITUENTS:
            raise ValueError("Unknown metal, backbone or phosphine substituent")

    @property
    def catalyst_id(self) -> str:
        return f"{self.metal}_{self.backbone}_{self.substituent}"


@dataclass
class GeneratedCatalyst:
    spec: CatalystSpec
    state: str
    atoms: Atoms
    donor_indices: tuple[int, int, int]
    proton_site_index: int
    ligand_indices: tuple[int, ...]
    ligand_smiles: str
    seed: int
    ligand_proton_index: int | None
    hydride_indices: tuple[int, ...]
    carbonyl_indices: tuple[tuple[int, int], ...]
    source_atom_maps: dict[int, int]
    embedding_status: str
    ligand_bond_indices: tuple[tuple[int, int], ...]
    metal_index: int = 0

    @property
    def charge(self) -> int:
        return 1 if self.state == "protonated_reference" else 0

    @property
    def unpaired_electrons(self) -> int:
        return 0

    def metadata(self) -> dict:
        symbols = self.atoms.get_chemical_symbols()
        return {
            **asdict(self.spec), "catalyst_id": self.spec.catalyst_id, "state": self.state,
            "formula": self.atoms.get_chemical_formula(), "charge": self.charge,
            "unpaired_electrons": self.unpaired_electrons,
            "spin_provenance": "Nominal low-spin occupation assumption; not a determined ground spin. GFN2-xTB is not a spin-state discriminating model.",
            "metal_oxidation_state": OXIDATION_STATES[self.spec.metal],
            "d_electrons": 8 if self.spec.metal == "Co" else 6,
            "nominal_valence_electrons": 18 if self.state == "hydrogenated" else 16,
            "metal_index": self.metal_index, "donor_indices": list(self.donor_indices),
            "donor_labels": [symbols[i] for i in self.donor_indices],
            "terminal_angle_label": "P-M-P" if self.spec.backbone == "macho_pnp" else "P-M-N",
            "P_M_P_applicable": self.spec.backbone == "macho_pnp",
            "proton_site_index": self.proton_site_index, "proton_site_element": symbols[self.proton_site_index],
            "proton_site_chemistry": {"macho_pnp": "amine/amido N-H", "pyridine_pnn": "phosphinomethyl C-H; aromatic carbanion resonance representation", "bipyridine_pnnoh": "hydroxypyridine O-H/O-; benzylic C-H is a second unscanned responsive site"}[self.spec.backbone],
            "ligand_proton_index": self.ligand_proton_index,
            "hydride_indices": list(self.hydride_indices),
            "carbonyl_indices": [list(pair) for pair in self.carbonyl_indices],
            "ligand_indices": list(self.ligand_indices),
            "ligand_bond_indices": [list(pair) for pair in self.ligand_bond_indices],
            "ligand_smiles_explicit_h": self.ligand_smiles,
            "atom_map_to_complex_index": self.source_atom_maps,
            "seed": self.seed, "embedding_status": self.embedding_status,
            "structure_status": "Designed starting geometry; requires electronic optimization and coordination-retention validation",
            "backbone_reference": BACKBONE_REFERENCES[self.spec.backbone],
            "literature_identity": "Scaffold-informed designed analogue; exact metal/R/ancillary combination is not asserted published",
            "transformations": "hydrogenated = active + H2 (one ligand proton and one metal hydride); protonated_reference = active + H+ (charge +1)",
            "energy_comparison_scope": "States of different stoichiometry/charge require explicit H2/proton/electron chemical-potential references",
            "coordinates_unit": "angstrom",
        }


def enumerate_specs() -> tuple[CatalystSpec, ...]:
    """Exactly 4 metals x 3 molecular backbones x 2 substituents = 24 designs."""
    return tuple(CatalystSpec(m, b, r) for m in METALS for b in BACKBONES for r in SUBSTITUENTS)


def _neutral_ligand(spec: CatalystSpec):
    from rdkit import Chem
    group = "c6ccccc6" if spec.substituent == "Ph" else "C(C)C"
    p1, p3 = f"[P:1]({group}){group}", f"[P:3]({group}){group}"
    if spec.backbone == "macho_pnp":
        smiles = f"[NH:2](CC{p1})CC{p3}"
        proton_map = 2
    elif spec.backbone == "pyridine_pnn":
        smiles = f"[CH2:4]({p1})c1cccc(C[N:3](CC)CC)[n:2]1"
        proton_map = 4
    else:
        smiles = f"[CH2:4]({p1})c1cccc(-c2cccc([OH:5])[n:3]2)[n:2]1"
        proton_map = 5
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"Internal ligand graph failed RDKit sanitization: {smiles}")
    mol = Chem.AddHs(mol)
    maps = {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum()}
    return mol, maps, proton_map


def _embed_ligand(mol, donors, seed: int):
    """Distance-geometry donor constraints followed by restrained UFF preparation.

    UFF is used for organic starting geometry only; its energy is neither saved
    nor interpreted as a quantum result. Metal atoms never enter UFF.
    """
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    targets = np.array([[2.30, -0.35, 0], [0, -2.05, 0], [-2.20, -0.35, 0]], dtype=float)
    for attempt in range(8):
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = (int(seed) + attempt * 7919) % (2**31 - 1)
        parameters.useRandomCoords = True
        parameters.SetCoordMap({i: Point3D(*p) for i, p in zip(donors, targets)})
        status = AllChem.EmbedMolecule(mol, parameters)
        if status != 0:
            continue
        coordinates = np.asarray(mol.GetConformer().GetPositions())
        source = coordinates[list(donors)]
        source_center, target_center = source.mean(axis=0), targets.mean(axis=0)
        u, _, vt = np.linalg.svd((source - source_center).T @ (targets - target_center))
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        coordinates = (coordinates - source_center) @ rotation + target_center
        conf = mol.GetConformer()
        for i, position in enumerate(coordinates):
            conf.SetAtomPosition(i, Point3D(*position))
        ff = AllChem.UFFGetMoleculeForceField(mol)
        if ff is None:
            raise RuntimeError("Organic ligand lacks UFF preparation parameters")
        for i, position in zip(donors, targets):
            extra = ff.AddExtraPoint(*position, fixed=True) - 1
            ff.AddDistanceConstraint(extra, i, 0, 0, 1500)
        # Reserve all ancillary sites across the metal/state matrix. These
        # geometric restraints prevent severe initial atom overlaps; they do
        # not constitute a metal force field or an electronic calculation.
        reserved = [(0, 0, -1.82), (0, 0, -2.96), (0, 0, 1.82),
                    (0, 0, 2.96), (0, 1.62, 0), (0, 0, 0)]
        for position in reserved:
            extra = ff.AddExtraPoint(*position, fixed=True) - 1
            for atom in mol.GetAtoms():
                lower = 1.25 if atom.GetAtomicNum() == 1 else 1.65
                ff.AddDistanceConstraint(extra, atom.GetIdx(), lower, 100.0, 300)
        ff.Initialize()
        ff.Minimize(maxIts=700)
        coordinates = np.asarray(conf.GetPositions())
        non_donor = [i for i in range(len(coordinates)) if i not in donors]
        if np.min(np.linalg.norm(coordinates[non_donor], axis=1)) > 1.20:
            return mol, f"ETKDG constrained embedding succeeded on attempt {attempt + 1}; organic UFF preparation only"
    raise RuntimeError("Constrained ligand embedding failed to avoid metal-center overlap after 8 seeds")


def generate_catalyst(spec: CatalystSpec, state: str = "active", *, seed: int = 20260913,
                      arm_torsion_deg: float = 0.0) -> GeneratedCatalyst:
    """Assemble an atom-resolved design with formal charge and atom mappings.

    ``arm_torsion_deg`` rotates a phosphine substituent around its P-C bond to
    seed a conformer search. Its value supplies a geometric seed, not a claim
    that a minimum/global minimum or a transition state has been located.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    state = {"active_deprotonated": "active", "active_amido": "active", "amido": "active", "hydride": "hydrogenated"}.get(state, state)
    if state not in STATES:
        raise ValueError(f"State must be one of {STATES}")
    if not np.isfinite(arm_torsion_deg):
        raise ValueError("arm_torsion_deg must be finite")
    mol, maps, proton_map = _neutral_ligand(spec)
    donors = tuple(maps[i] for i in (1, 2, 3))
    mol, embedding_status = _embed_ligand(mol, donors, seed)
    if arm_torsion_deg:
        p = mol.GetAtomWithIdx(maps[1])
        choices = [atom for atom in p.GetNeighbors() if atom.GetAtomicNum() == 6]
        for carbon in choices:
            outer = [a for a in carbon.GetNeighbors() if a.GetIdx() != p.GetIdx()]
            inner = [a for a in p.GetNeighbors() if a.GetIdx() != carbon.GetIdx()]
            path_to_central_donor = Chem.GetShortestPath(mol, carbon.GetIdx(), donors[1])
            if outer and inner and p.GetIdx() in path_to_central_donor:
                indexes = (inner[0].GetIdx(), p.GetIdx(), carbon.GetIdx(), outer[0].GetIdx())
                current = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *indexes)
                rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *indexes, current + arm_torsion_deg)
                break
    site = maps[proton_map]
    attached_h = [atom.GetIdx() for atom in mol.GetAtomWithIdx(site).GetNeighbors() if atom.GetAtomicNum() == 1]
    if not attached_h:
        raise RuntimeError("Defined proton-responsive site has no removable hydrogen")
    proton = attached_h[0]
    if state == "active":
        editable = Chem.RWMol(mol)
        editable.RemoveAtom(proton)
        atom = editable.GetAtomWithIdx(site)
        atom.SetFormalCharge(-1)
        atom.SetNoImplicit(True)
        mol = editable.GetMol()
        Chem.SanitizeMol(mol)
        proton = None
    positions = [np.zeros(3)] + [p for p in mol.GetConformer().GetPositions()]
    symbols = [spec.metal] + [atom.GetSymbol() for atom in mol.GetAtoms()]
    ligand_indices = tuple(range(1, len(symbols)))
    carbonyls, hydrides = [], []
    # Nominal 16e sets: five-coordinate Ru/Mn/Fe and four-coordinate Co.
    co_directions = [(0, 0, -1)]
    if spec.metal == "Mn":
        co_directions.append((0, 0, 1))
    for direction in co_directions:
        vector = np.array(direction, dtype=float)
        carbonyls.append((len(symbols), len(symbols) + 1))
        symbols.extend(["C", "O"])
        positions.extend([1.82 * vector, 2.96 * vector])
    if spec.metal in ("Ru", "Fe"):
        hydrides.append(len(symbols))
        symbols.append("H")
        positions.append(np.array([0, 0, 1.65]))
    if state == "hydrogenated":
        hydrides.append(len(symbols))
        symbols.append("H")
        positions.append(np.array([0, 1.62, 0]))
    atoms = Atoms(symbols=symbols, positions=positions)
    distances = atoms.get_all_distances() + np.eye(len(atoms)) * 1e6
    if np.min(distances) < 0.75:
        raise ValueError("Generated conformer has an interatomic separation below 0.75 angstrom; try another seed or torsion")
    result = GeneratedCatalyst(spec=spec, state=state, atoms=atoms,
        donor_indices=tuple(i + 1 for i in donors), proton_site_index=site + 1,
        ligand_indices=ligand_indices, ligand_smiles=Chem.MolToSmiles(mol), seed=seed,
        ligand_proton_index=None if proton is None else proton + 1,
        hydride_indices=tuple(hydrides), carbonyl_indices=tuple(carbonyls),
        source_atom_maps={key: value + 1 for key, value in maps.items()}, embedding_status=embedding_status,
        ligand_bond_indices=tuple((bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1) for bond in mol.GetBonds()))
    if (sum(atoms.numbers) - result.charge - result.unpaired_electrons) % 2:
        raise ValueError("Electronic occupation parity is inconsistent with generated formula and charge")
    atoms.info.update(result.metadata())
    return result


def write_library(output_dir: str | Path, *, states=STATES, seed: int = 20260913) -> list[dict]:
    """Write XYZ structures plus JSON identities; no calculated energies are emitted."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for index, spec in enumerate(enumerate_specs()):
        for state in states:
            generated = generate_catalyst(spec, state, seed=seed + index)
            filename = f"{spec.catalyst_id}__{generated.state}.xyz"
            write(output_dir / filename, generated.atoms, format="xyz")
            record = {**generated.metadata(), "xyz_file": filename}
            (output_dir / filename.replace(".xyz", ".json")).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
            records.append(record)
    (output_dir / "library_index.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    return records


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--states", nargs="+", choices=STATES, default=list(STATES))
    parser.add_argument("--seed", type=int, default=20260913)
    args = parser.parse_args()
    records = write_library(args.output_dir, states=args.states, seed=args.seed)
    print(f"Generated {len(records)} designed starting structures; no quantum energies computed")
