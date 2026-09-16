"""Read-only-source import of the published Mn PNP complex-1 crystal structure.

Uses gemmi for CIF/symmetry and ASE for cell conversion/periodic contacts.
Extracted structures remain unoptimized crystal candidates, not solution labels.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil

import gemmi
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.geometry import cellpar_to_cell
from ase.io import write
from ase.neighborlist import neighbor_list

EXPECTED_SHA256 = "9704aa1897177a1026ea7a87f57cf2679bd0bd5c6647910ef60cc06696b663cd"
EXPECTED_FORMULA = Counter(C=18, H=37, Br=1, Mn=1, N=1, O=2, P=2)
SOURCE_URL = "https://ndownloader.figshare.com/files/5467289"
SOURCE_DOI = "10.1021/jacs.6b03709.s002"
SOURCE_API = "https://api.figshare.com/v2/articles/3472229"
DATA_LICENSE = {"name": "CC BY-NC 4.0", "url": "https://creativecommons.org/licenses/by-nc/4.0/",
                "verified_from": SOURCE_API, "verified_date": "2026-09-16"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def normalized_cif(text):
    """Repair only this known preamble, preserving all crystallographic fields."""
    if re.search(r"(?m)^data_", text):
        return text, {"changed": False}
    position = text.find("_audit_creation_method")
    if position < 0 or text[:position].strip() != "CIF files of complex 1":
        raise ValueError("Unsupported CIF preamble; refusing to guess a data block")
    return "data_published_complex_1\n" + text[position:], {
        "changed": True, "operation": "Replace only non-CIF introductory text with a data block identifier",
        "removed_preamble": text[:position], "added_header": "data_published_complex_1"}


def formula_counts(text):
    compact = gemmi.cif.as_string(text).replace(" ", "")
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", compact)
    if "".join(a + n for a, n in parts) != compact:
        raise ValueError("Unsupported formula syntax")
    counts = Counter()
    for symbol, number in parts:
        counts[symbol] += int(number or 1)
    return counts


def parse_structure(text):
    block = gemmi.cif.read_string(text).sole_block()
    tags = ["label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
            "disorder_assembly", "disorder_group", "calc_flag"]
    table = block.find("_atom_site_", tags)
    if not len(table):
        raise ValueError("Missing full site/occupancy/disorder loop")
    sites = []
    for row in table:
        site = dict(label=row[0], symbol=row[1], fractional=[gemmi.cif.as_number(row[i]) for i in (2, 3, 4)],
                    occupancy=gemmi.cif.as_number(row[5]), disorder_assembly=row[6], disorder_group=row[7],
                    calculation_flag=row[8])
        if site["symbol"] not in atomic_numbers or not np.all(np.isfinite(site["fractional"] + [site["occupancy"]])):
            raise ValueError("Nonfinite or unknown atom site")
        sites.append(site)
    if len({site["label"] for site in sites}) != len(sites):
        raise ValueError("Duplicate atom labels")
    cellpar = [gemmi.cif.as_number(block.find_value("_cell_" + key)) for key in
               ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")]
    if not np.all(np.isfinite(cellpar)) or min(cellpar[:3]) <= 0:
        raise ValueError("Invalid unit cell")
    cell = cellpar_to_cell(cellpar)
    if np.linalg.det(cell) <= 0:
        raise ValueError("Non-positive unit cell")
    operations = [gemmi.cif.as_string(value) for value in block.find_values("_space_group_symop_operation_xyz")]
    if not operations:
        raise ValueError("Explicit CIF symmetry operations required")
    bonds = []
    for row in block.find("_geom_bond_", ["atom_site_label_1", "atom_site_label_2", "distance", "site_symmetry_2"]):
        if row[3] not in (".", "1_555"):
            raise ValueError("This importer requires identity-symmetry published bonds; do not guess a periodic bond")
        bonds.append((row[0], row[1], gemmi.cif.as_number(row[2])))
    if not bonds:
        raise ValueError("Published bond table required")
    return dict(sites=sites, cell=cell, cell_parameters=cellpar, symmetry_operations=operations, bonds=bonds,
                formula=formula_counts(block.find_value("_chemical_formula_sum")),
                Z=int(gemmi.cif.as_number(block.find_value("_cell_formula_units_Z"))),
                temperature_K=gemmi.cif.as_number(block.find_value("_cell_measurement_temperature")),
                hydrogen_treatment=gemmi.cif.as_string(block.find_value("_refine_ls_hydrogen_treatment")),
                R1_all=gemmi.cif.as_number(block.find_value("_refine_ls_R_factor_all")),
                R1_observed=gemmi.cif.as_number(block.find_value("_refine_ls_R_factor_gt")))


def expand_cell(parsed):
    provenance, fractions, symbols = [], [], []
    for operation_index, text in enumerate(parsed["symmetry_operations"]):
        operation = gemmi.Op(text)
        for site_index, site in enumerate(parsed["sites"]):
            fractional = np.mod(operation.apply_to_xyz(site["fractional"]), 1.)
            for old_symbol, old_fractional in zip(symbols, fractions):
                delta = fractional - old_fractional
                delta -= np.round(delta)
                if np.linalg.norm(delta @ parsed["cell"]) < 1e-5:
                    raise ValueError("Special-position/overlapping symmetry sites require dedicated occupancy handling")
            fractions.append(fractional)
            symbols.append(site["symbol"])
            provenance.append(dict(site_index=site_index, label=site["label"], symmetry_operation_index=operation_index))
    atoms = Atoms(symbols=symbols, scaled_positions=fractions, cell=parsed["cell"], pbc=True)
    return atoms, provenance


def periodic_components(atoms, scale=1.20):
    """Build/unwrap finite molecular components; reject nonzero lattice cycles."""
    ii, jj, shifts = neighbor_list("ijS", atoms, [scale * covalent_radii[z] for z in atoms.numbers], self_interaction=False)
    graph = [[] for _ in atoms]
    for i, j, shift in zip(ii, jj, shifts):
        graph[int(i)].append((int(j), np.asarray(shift, dtype=int)))
    components, seen, infinite = [], set(), []
    for first in range(len(atoms)):
        if first in seen:
            continue
        offsets = {first: np.zeros(3, dtype=int)}
        queue = deque([first])
        while queue:
            i = queue.popleft()
            seen.add(i)
            for j, shift in graph[i]:
                proposed = offsets[i] + shift
                if j in offsets:
                    if not np.array_equal(offsets[j], proposed):
                        infinite.append(dict(i=i, j=j, residual=(proposed - offsets[j]).tolist()))
                else:
                    offsets[j] = proposed
                    queue.append(j)
        indices = sorted(offsets)
        components.append(dict(indices=indices, offsets=[offsets[i].tolist() for i in indices]))
    return components, graph, infinite


def check_published_bonds(parsed):
    sites = {s["label"]: s for s in parsed["sites"]}
    rows = []
    for a, b, distance in parsed["bonds"]:
        if a not in sites or b not in sites or not np.isfinite(distance):
            raise ValueError("Invalid published bond")
        difference = np.asarray(sites[a]["fractional"]) - sites[b]["fractional"]
        calculated = float(np.linalg.norm(difference @ parsed["cell"]))
        rows.append(dict(atom_1=a, atom_2=b, reported_distance_angstrom=distance,
                         recalculated_distance_angstrom=calculated,
                         absolute_deviation_angstrom=abs(calculated - distance)))
    return rows


def molecule_metadata(atoms, provenance, edges, molecule_id):
    symbols = atoms.get_chemical_symbols()
    metal = symbols.index("Mn")
    neighbors = {i: set() for i in range(len(atoms))}
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    donors = sorted(i for i in neighbors[metal] if symbols[i] in {"P", "N"})
    nitrogens = [i for i in donors if symbols[i] == "N"]
    if len(nitrogens) != 1:
        raise ValueError("Exactly one coordinated N required")
    n = nitrogens[0]
    nh = [i for i in neighbors[n] if symbols[i] == "H"]
    if len(nh) != 1:
        raise ValueError("Exactly one explicit NH hydrogen required")
    carbonyls = []
    for c in neighbors[metal]:
        if symbols[c] == "C":
            oxygen = [i for i in neighbors[c] if symbols[i] == "O"]
            if len(oxygen) != 1:
                raise ValueError("Metal-bound carbon is not a unique CO")
            carbonyls.append([c, oxygen[0]])
    halides = [i for i in neighbors[metal] if symbols[i] == "Br"]
    if Counter(symbols[i] for i in neighbors[metal]) != Counter(P=2, N=1, C=2, Br=1):
        raise ValueError("Expected octahedral P2 N C2 Br coordination not found")
    excluded = {metal, *halides, *(i for pair in carbonyls for i in pair)}
    ligand = [i for i in range(len(atoms)) if i not in excluded]
    return dict(catalyst_id="Mn_PNP_NH_iPr_complex1_" + molecule_id, metal="Mn", backbone="literature_PNP_NH",
                substituent="iPr", state="precatalyst", formula=atoms.get_chemical_formula(),
                charge=0, multiplicity=1, unpaired_electrons=0,
                charge_provenance="Neutral molecular precursor working assignment from neutral PNP/CO and bound bromide; review before quantum submission",
                spin_provenance="Singlet is a chemical starting candidate, not an experimental spin-state measurement",
                metal_oxidation_state=1, oxidation_state_provenance="Formal electron-counting hypothesis; not independently measured by this import",
                metal_index=metal, donor_indices=donors, donor_labels=[symbols[i] for i in donors],
                proton_site_index=n, proton_site_element="N", proton_site_chemistry="Secondary-amine NH in published precursor",
                ligand_proton_index=nh[0], ligand_indices=ligand,
                ligand_bond_indices=[[a, b] for a, b in edges if a in ligand and b in ligand],
                carbonyl_indices=carbonyls, halide_indices=halides, hydride_indices=[],
                all_bond_indices=edges, source_atom_provenance=provenance,
                coordinate_provenance="Unrelaxed single-crystal X-ray geometry at 150 K; original H positions retained",
                source_doi=SOURCE_DOI, source_cif_sha256=EXPECTED_SHA256,
                structure_origin="public_crystal_derived_candidate",
                provenance={"source_url": SOURCE_URL, "source_identifier": "doi:" + SOURCE_DOI,
                            "structure_derivation": "Preserved original CIF; repaired only missing data header in a copy; expanded all symmetry operations; unwrapped finite periodic molecular graph; selected one of two independent crystal molecules; translated Mn to origin; retained all original H sites; no geometry optimization"},
                third_party_data_license=DATA_LICENSE,
                active_state_generated=False, geometry_optimized=False, chemistry_review_required=True,
                compatible_with_generated_catalog_identity_function=False)


def import_structure(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Existing import is immutable; supply a fresh output directory")
    if digest(source) != EXPECTED_SHA256:
        raise ValueError("Source CIF checksum differs from the published download; refuse substituted structure")
    raw = source.read_text(encoding="utf-8")
    normalized, repair = normalized_cif(raw)
    parsed = parse_structure(normalized)
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(source, output / "original.cif")
    (output / "parsed_copy.cif").write_text(normalized, encoding="utf-8")
    dump(output / "normalization.json", repair)
    shutil.copyfile(Path(__file__), output / "importer_at_execution.py")
    report = dict(schema="public_catalyst_import_v1", started_utc=datetime.now(timezone.utc).isoformat(),
                  status="blocked", source_url=SOURCE_URL, source_doi=SOURCE_DOI, source_sha256=digest(source),
                  official_metadata_api=SOURCE_API, third_party_data_license=DATA_LICENSE,
                  parser="gemmi " + gemmi.__version__, original_modified=False,
                  normalized_sha256=digest(output / "parsed_copy.cif"), blockers=[], molecules=[],
                  crystallographic_formula=dict(parsed["formula"]), cell_parameters=parsed["cell_parameters"],
                  Z=parsed["Z"], asymmetric_unit_atom_count=len(parsed["sites"]),
                  temperature_K=parsed["temperature_K"], hydrogen_treatment=parsed["hydrogen_treatment"],
                  R1_all=parsed["R1_all"], R1_observed=parsed["R1_observed"])
    try:
        partial = [s["label"] for s in parsed["sites"] if abs(s["occupancy"] - 1.) > 1e-8]
        disorder = [s["label"] for s in parsed["sites"] if s["disorder_assembly"] not in {".", "?"} or s["disorder_group"] not in {".", "?"}]
        report.update(partial_occupancy_sites=partial, disorder_sites=disorder,
                      atom_site_calculation_flags=dict(Counter(s["calculation_flag"] for s in parsed["sites"])))
        if partial or disorder:
            raise ValueError("Partial occupancy/disorder requires explicit alternate-conformer handling")
        if parsed["formula"] != EXPECTED_FORMULA:
            raise ValueError("Published formula differs from target complex 1")
        atoms, provenance = expand_cell(parsed)
        components, graph, infinite = periodic_components(atoms)
        report.update(unit_cell_atom_count=len(atoms), symmetry_operation_count=len(parsed["symmetry_operations"]),
                      periodic_component_count=len(components), infinite_periodic_cycles=infinite,
                      connectivity_rule="1.20 times sum of ASE covalent radii, periodic minimum-image contacts",
                      components=[dict(atom_count=len(c["indices"]), formula=dict(Counter(atoms[i].symbol for i in c["indices"])),
                                       crosses_cell_boundary=any(any(v) for v in c["offsets"])) for c in components])
        if infinite or len(components) != parsed["Z"] or any(Counter(atoms[i].symbol for i in c["indices"]) != EXPECTED_FORMULA for c in components):
            raise ValueError("Full symmetry-expanded cell does not yield Z complete finite target molecules")
        bond_check = check_published_bonds(parsed)
        dump(output / "published_bond_checks.json", bond_check)
        report["maximum_bond_distance_reconstruction_error_angstrom"] = max(row["absolute_deviation_angstrom"] for row in bond_check)
        if report["maximum_bond_distance_reconstruction_error_angstrom"] > .002:
            raise ValueError("Reconstructed coordinates do not reproduce published bond distances within rounding tolerance")
        published_edges = {tuple(sorted((a, b))) for a, b, distance in parsed["bonds"]}
        selected = []
        for component in components:
            metal_indices = [i for i in component["indices"] if atoms[i].symbol == "Mn"]
            if len(metal_indices) != 1:
                raise ValueError("Nonmononuclear component")
            if provenance[metal_indices[0]]["symmetry_operation_index"] == 0:
                selected.append((component, metal_indices[0]))
        if len(selected) != 2:
            raise ValueError("Expected two crystallographically independent precursor molecules")
        for component, metal_global in selected:
            global_indices = component["indices"]
            offsets = {i: np.array(v) for i, v in zip(global_indices, component["offsets"])}
            ordered = [metal_global] + [i for i in global_indices if i != metal_global]
            lookup = {i: j for j, i in enumerate(ordered)}
            positions = np.array([atoms.positions[i] + offsets[i] @ atoms.cell.array for i in ordered])
            positions -= positions[0]
            molecule = Atoms([atoms[i].symbol for i in ordered], positions=positions)
            edges = sorted({tuple(sorted((lookup[i], lookup[j]))) for i in ordered for j, shift in graph[i] if j in lookup})
            actual_edges = {tuple(sorted((provenance[ordered[a]]["label"], provenance[ordered[b]]["label"]))) for a, b in edges}
            labels = {provenance[i]["label"] for i in ordered}
            expected_edges = {edge for edge in published_edges if set(edge) <= labels}
            if actual_edges != expected_edges:
                raise ValueError("Distance-derived molecular graph differs from published CIF bond list")
            molecule_id = provenance[metal_global]["label"]
            destination = output / molecule_id
            destination.mkdir()
            metadata = molecule_metadata(molecule, [provenance[i] for i in ordered], [list(e) for e in edges], molecule_id)
            write(destination / "crystal_candidate.xyz", molecule, format="xyz")
            metadata["geometry_sha256"] = digest(destination / "crystal_candidate.xyz")
            dump(destination / "metadata.json", metadata)
            metal_bonds = [dict(index=i, symbol=molecule[i].symbol, source_label=provenance[ordered[i]]["label"],
                               distance_angstrom=float(molecule.get_distance(0, i))) for i in range(1, len(molecule))
                           if (0, i) in edges]
            report["molecules"].append(dict(id=molecule_id, atoms=len(molecule), formula=molecule.get_chemical_formula(),
                xyz=(destination / "crystal_candidate.xyz").relative_to(output).as_posix(),
                xyz_sha256=metadata["geometry_sha256"], metadata_sha256=digest(destination / "metadata.json"),
                NH_distance_angstrom=float(molecule.get_distance(metadata["proton_site_index"], metadata["ligand_proton_index"])),
                metal_neighbors=metal_bonds, graph_matches_published_bonds=True))
        report.update(status="structure_extraction_verified", complete_molecule_verified=True,
                      quantum_submission_ready=False, solution_catalyst_identity_validated=False,
                      charge_spin_experimentally_verified=False,
                      next_gate="Chemist reviews precursor identity, charge/spin candidates and explicit activation chemistry before DFT submission")
    except Exception as error:
        report["blockers"].append(f"{type(error).__name__}: {error}")
    finally:
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        dump(output / "structure_verification.json", report)
        dump(output / "file_manifest.json", {str(p.relative_to(output)).replace("\\", "/"): digest(p)
             for p in sorted(output.rglob("*")) if p.is_file() and p.name != "file_manifest.json"})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "data/phase4/public_structure")
    args = parser.parse_args()
    report = import_structure(args.source, args.output)
    print(json.dumps(report, indent=2))
    if report["status"] != "structure_extraction_verified":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
