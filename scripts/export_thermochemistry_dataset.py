"""Export only source-verified, identity-retaining actual minimum certificates.

Run from any checkout: historical absolute paths remain in original JSON, while
local same-directory evidence is preferred and all exports use repository paths.
The chosen geometry minimizes certified G at 383.15 K, not the electronic-energy
search index; its own full 13-temperature series is used without geometry swaps.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import gzip
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import sys

os.environ.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import numpy as np
from ase import Atoms, units
from ase.io import read
from pincer_catmech.generators.combinatorial_pincer import enumerate_specs, STATES
from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.thermochemistry import DEFAULT_TEMPERATURES, project_hessian, thermochemistry_from_modes
from pincer_catmech.quantum.xtb_backend import atomic_json, parse_energy, utc_now, validate_electrons

BASE_TEMPERATURE = 383.15
CYCLE_TOLERANCE_KCAL_MOL = 1e-8
ROOTS = ("data/campaign/thermochemistry", "data/campaign/thermochemistry_identity_replacements",
         "data/targeted_rescue/thermochemistry")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return Path(path).resolve().relative_to(REPO).as_posix()


def historical_key(path):
    return str(path).replace("\\", "/").rstrip("/").casefold()


def resolve_evidence(certificate, historical_source, expected_sha256):
    """Prefer archived local evidence; a filename alone never establishes identity."""
    if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("missing or malformed evidence SHA256")
    name = PurePosixPath(str(historical_source).replace("\\", "/")).name
    candidates = [certificate.parent / name, REPO / str(historical_source).replace("\\", "/"), Path(historical_source)]
    for candidate in candidates:
        if candidate.is_file() and digest(candidate) == expected_sha256:
            candidate.resolve().relative_to(REPO)
            return candidate.resolve()
    raise ValueError(f"no local evidence with matching SHA256: {name}")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference_identity(atoms, metadata, charge):
    from rdkit import Chem
    from run_reference_thermochemistry import connectivity_check
    molecule = Chem.AddHs(Chem.MolFromSmiles(metadata["smiles"]))
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    formal_charge = Chem.GetFormalCharge(molecule)
    if metadata.get("state_model") == "contact ion pair":
        symbols.append("K")
        formal_charge += 1
    require(symbols == atoms.get_chemical_symbols() and formal_charge == charge, "reference atom order/formal charge does not match its molecular graph")
    return connectivity_check(atoms, [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in molecule.GetBonds()])


@dataclass
class CertifiedMinimum:
    key: str
    kind: str
    certificate: Path
    geometry: Path
    spectrum: Path
    atoms: Atoms
    charge: int
    data: dict
    metadata: dict
    rows: dict
    ledger: dict


def validate_certificate(path: Path, history: dict, *, reference=False):
    raw = path.read_bytes()
    try:
        data = json.loads(raw)
        require(isinstance(data, dict), "certificate must be a JSON object")
    except (ValueError, TypeError) as error:
        return {"key": path.parent.name.replace("__", "|"), "kind": "reference" if reference else "catalyst",
                "certificate": relative(path), "certificate_sha256": hashlib.sha256(raw).hexdigest(),
                "reported_status": None, "eligible": False, "validation_status": "malformed_certificate",
                "reason": f"{type(error).__name__}: {error}"}, None
    source = data.get("source_metadata", {})
    key = (source.get("index_key") or (f"{source.get('catalyst_id')}|{source.get('state')}" if source.get("catalyst_id")
           else source.get("species", path.parent.name))) if isinstance(source, dict) else path.parent.name.replace("__", "|")
    ledger = {"key": key, "kind": "reference" if reference else "catalyst", "certificate": relative(path),
              "certificate_sha256": hashlib.sha256(raw).hexdigest(), "reported_status": data.get("status"),
              "eligible": False, "validation_status": "not_checked"}
    try:
        require(isinstance(source, dict), "source_metadata must be a JSON object")
        require(data.get("finished_utc"), "calculation_in_progress")
        require(data.get("accepted") is True and data.get("status") == "verified_minimum", "full minimum was not certified")
        require(data.get("method") == "GFN2-xTB" and data.get("solvent") == "toluene" and data.get("solvation_state") == "gsolv", "electronic method or solvent convention differs")
        require(data.get("electronic_temperature_K") == 300 and not data.get("transition_state"), "unexpected electronic temperature or TS certificate")
        geometry = path.parent / "accepted_geometry.xyz"
        require(digest(geometry) == data["accepted_geometry_sha256"], "accepted geometry SHA256 mismatch")
        atoms = read(geometry)
        charge, unpaired = data["charge"], data["unpaired_electrons"]
        validate_electrons(atoms, charge, unpaired)
        accepted = data["accepted_attempt"]
        require(isinstance(accepted, int) and not isinstance(accepted, bool) and accepted >= 0, "invalid accepted attempt index")
        attempt = data["attempts"][accepted]
        require(attempt["status"] == "verified_minimum", "accepted attempt is not a minimum")
        require(0 <= attempt["max_force_eV_A"] <= min(.03, data["force_threshold_eV_A"]), "true-force convergence gate failed")
        require(0 <= attempt["hessian_antisymmetry_relative"] <= .05 and 0 <= attempt["gradient_hessian_energy_difference_eV"] <= .01, "Hessian consistency gate failed")
        records = data["thermochemistry"]
        temperatures = [row["temperature"] for row in records]
        require(len(temperatures) == len(DEFAULT_TEMPERATURES) and set(temperatures) == set(DEFAULT_TEMPERATURES), "missing or duplicate temperature rows")
        spectrum = resolve_evidence(path, records[0]["source"], records[0]["source_sha256"])
        npz_sha = digest(spectrum)
        require(all(row["source_sha256"] == npz_sha and PurePosixPath(row["source"].replace("\\", "/")).name == spectrum.name for row in records), "temperature rows use different Hessian sources")
        with np.load(spectrum, allow_pickle=False) as npz:
            require(np.array_equal(atoms.numbers, npz["atomic_numbers"]), "Hessian atom mapping differs from accepted geometry")
            require(np.allclose(atoms.positions, npz["positions"], atol=1e-7, rtol=0), "Hessian coordinates differ from accepted geometry")
            atoms.positions = npz["positions"]
            atoms.set_masses(npz["masses"])
            matrix = npz["hessian_ev_a2"].copy()
            projected = project_hessian(atoms, matrix)
            frequencies = npz["frequencies_cm1"].copy()
            require(projected.external_rank == int(npz["external_rank"]) == attempt["external_modes_projected"], "external-mode projection ranks disagree")
        require(np.all(np.isfinite(frequencies)) and np.all(frequencies > 0), "nonpositive or nonfinite internal modes")
        require(len(frequencies) == 3 * len(atoms) - projected.external_rank, "not a full internal spectrum")
        require(np.allclose(frequencies, projected.frequencies_cm1, atol=1e-4, rtol=1e-7) and
                np.allclose(frequencies, attempt["frequencies_cm1"], atol=1e-7, rtol=0) and attempt["imaginary_count"] == 0,
                "stored frequencies disagree with the actual full Hessian")
        trace = path.parent / f"gradient_{accepted:02d}.jsonl.gz"
        gradients = [json.loads(line) for line in gzip.decompress(trace.read_bytes()).decode("utf-8").splitlines() if line.strip()]
        force_record = gradients[-1]
        forces = np.asarray(force_record["forces_eV_A"], dtype=float)
        require(forces.shape == (len(atoms), 3) and np.all(np.isfinite(forces)), "invalid actual gradient trace")
        require(force_record["symbols"] == atoms.get_chemical_symbols() and np.allclose(force_record["positions_A"], atoms.positions, atol=1e-7, rtol=0), "gradient trace belongs to a different geometry")
        require(abs(float(np.linalg.norm(forces, axis=1).max()) - attempt["max_force_eV_A"]) < 1e-8, "true-force trace and convergence certificate disagree")
        require(force_record["solvent"] == "toluene" and force_record["solvation_state"] == "gsolv" and force_record["charge"] == charge, "gradient electronic state differs")
        require(force_record["unpaired_electrons"] == unpaired and force_record["electronic_temperature_K"] == 300 and
                abs(force_record["energy_eV"] - attempt["gradient_energy_eV"]) < 1e-8, "gradient occupation/energy record differs")
        argv = attempt["native_command"]["argv"]
        for flag, value in (("--gfn", "2"), ("--chrg", str(charge)), ("--uhf", str(unpaired))):
            require(argv.count(flag) == 1 and argv[argv.index(flag)+1] == value, f"native command differs: {flag}")
        require(argv.count("--alpb") == 1 and argv[argv.index("--alpb")+1:argv.index("--alpb")+3] == ["toluene", "gsolv"] and "--hess" in argv, "native Hessian solvent command differs")
        raw_hessian = attempt["unprojected_hessian"]
        raw_path = resolve_evidence(path, raw_hessian["file"], raw_hessian["gzip_sha256"])
        raw_bytes = gzip.decompress(raw_path.read_bytes())
        require(hashlib.sha256(raw_bytes).hexdigest() == raw_hessian["source_sha256"], "raw Cartesian Hessian SHA256 mismatch")
        raw_matrix = np.array([float(token.replace("D", "E").replace("d", "e")) for token in raw_bytes.decode("utf-8").split()])
        require(raw_matrix.size == (3*len(atoms))**2 and np.all(np.isfinite(raw_matrix)), "raw native Hessian cardinality/finiteness failed")
        raw_matrix = raw_matrix.reshape(matrix.shape) * units.Hartree / units.Bohr**2
        require(np.allclose(matrix, (raw_matrix + raw_matrix.T)/2, atol=1e-8, rtol=1e-8), "NPZ Hessian differs from actual native Cartesian matrix")
        native = attempt["native_output"]
        native_path = resolve_evidence(path, native["file"], native["gzip_sha256"])
        native_bytes = gzip.decompress(native_path.read_bytes())
        require(hashlib.sha256(native_bytes).hexdigest() == native["source_sha256"], "native output content SHA256 mismatch")
        text = native_bytes.decode("utf-8", errors="replace")
        require("normal termination of xtb" in text.lower() and "Numerical Hessian" in text and
                re.search(r"frozen atoms in %\s*:\s*0\.0+\s+0", text) and
                re.search(r"Hessian scale factor\s*:\s*1\.0+", text), "native full/unscaled/unfrozen Hessian evidence missing")
        require(abs(parse_energy(text) - attempt["energy_eV"]) < 1e-8 and
                abs(force_record["energy_eV"] - attempt["energy_eV"]) <= .01, "native and recorded electronic energies disagree")
        if reference:
            metadata = source
            identity = reference_identity(atoms, metadata, charge)
        else:
            metadata = data.get("chemical_identity_metadata") or source.get("metadata")
            if metadata is None:
                source_directory = PurePosixPath(source["native_geometry"].replace("\\", "/")).parent
                metadata = history[historical_key(source_directory)]
            key = f"{metadata['catalyst_id']}|{metadata['state']}"
            ledger["key"] = key
            require(metadata["charge"] == charge and metadata["unpaired_electrons"] == unpaired, "identity electronic-state metadata differs")
            identity = catalyst_identity_screen(atoms, metadata)
            donor_distances = atoms.get_distances(metadata["metal_index"], metadata["donor_indices"])
            require(np.all((donor_distances > 1.4) & (donor_distances < 3.1)), "pincer donor coordination screen failed")
            ledger["metal_donor_distances_A"] = donor_distances.tolist()
        ledger["chemical_identity"] = identity
        require(identity["retained"], "intended molecular identity is not retained")
        recalculated = thermochemistry_from_modes(atoms, attempt["energy_eV"], frequencies,
            unpaired=unpaired, symmetry_number=data["symmetry_number"], temperatures=DEFAULT_TEMPERATURES)
        by_temperature = {row["temperature"]: row for row in records}
        for expected in recalculated:
            row = by_temperature[expected.temperature]
            require(row["energy_unit"] == "kcal/mol" and row["entropy_unit"] == "cal/(mol K)" and
                    row["pressure_pa"] == 101325 and row["concentration_mol_l"] == 1 and row["cutoff_cm1"] == 100 and
                    row["rotor_inertia"] == "grimme" and row["program"] == "GFN2-xTB+ALPB", "thermochemistry protocol or units differ")
            require(len(row["frequencies_cm1"]) == len(frequencies) and np.allclose(row["frequencies_cm1"], frequencies, atol=1e-7, rtol=0), "temperature row mode list differs from source spectrum")
            for field in ("E_elec", "ZPVE", "H_298", "S_harmonic", "S_qRRHO", "G_298_qRRHO_sol", "delta_G_concentration"):
                require(math.isfinite(row[field]) and abs(row[field] - getattr(expected, field)) < 1e-4,
                        f"thermochemistry no longer reproduces from source: {field} at {expected.temperature} K")
        ledger.update(eligible=True, validation_status="verified", accepted_geometry=relative(geometry),
            accepted_geometry_sha256=digest(geometry), spectrum=relative(spectrum), spectrum_sha256=npz_sha,
            gradient_trace=relative(trace), gradient_trace_sha256=digest(trace), native_output=relative(native_path),
            native_output_sha256=digest(native_path), G_383p15_kcal_mol=by_temperature[BASE_TEMPERATURE]["G_298_qRRHO_sol"],
            minimum_frequency_cm1=float(frequencies.min()), internal_mode_count=len(frequencies),
            accepted_attempt=accepted, max_force_eV_A=attempt["max_force_eV_A"])
        return ledger, CertifiedMinimum(key, ledger["kind"], path, geometry, spectrum, atoms, charge, data, metadata, by_temperature, ledger)
    except Exception as error:
        ledger.update(validation_status="calculation_in_progress" if not data.get("finished_utc") else "rejected", reason=f"{type(error).__name__}: {error}")
        return ledger, None


def select_certificates(certificates):
    selected = {}
    for certificate in certificates:
        previous = selected.get(certificate.key)
        if previous is None or (certificate.rows[BASE_TEMPERATURE]["G_298_qRRHO_sol"], relative(certificate.certificate)) < (previous.rows[BASE_TEMPERATURE]["G_298_qRRHO_sol"], relative(previous.certificate)):
            selected[certificate.key] = certificate
    return selected


def reaction_balance(reactants, products):
    counts = Counter()
    charge = 0
    for sign, side in ((-1, reactants), (1, products)):
        for species in side:
            counts.update({element: sign * number for element, number in Counter(species.atoms.get_chemical_symbols()).items()})
            charge += sign * species.charge
    difference = {element: number for element, number in sorted(counts.items()) if number}
    return {"balanced": not difference and charge == 0, "product_minus_reactant_atoms": difference, "product_minus_reactant_charge": charge}


def reaction_delta(reactants, products, temperature):
    return math.fsum([c.rows[temperature]["G_298_qRRHO_sol"] for c in products] +
                     [-c.rows[temperature]["G_298_qRRHO_sol"] for c in reactants])


def species_stoichiometry(steps):
    """Keep signed coefficients; equal molecular formula is not species identity."""
    total = Counter()
    for lhs, rhs in steps:
        total.update(rhs)
        total.subtract(lhs)
    return {key: count for key, count in sorted(total.items()) if count}


def export_cycle(selected, identifier):
    active, hydrogenated = f"{identifier}|active", f"{identifier}|hydrogenated"
    steps = {
        "alcohol_dehydrogenation": ([active, "benzyl_alcohol"], [hydrogenated, "benzaldehyde"]),
        "hemiaminal_formation": (["benzaldehyde", "aniline"], ["hemiaminal"]),
        "dehydration": (["hemiaminal"], ["imine", "water"]),
        "imine_hydrogenation": ([hydrogenated, "imine"], [active, "product_amine"]),
    }
    net = (["benzyl_alcohol", "aniline"], ["product_amine", "water"])
    keys = sorted({key for pair in steps.values() for side in pair for key in side})
    status = {"catalyst_id": identifier, "steps": steps, "net_reaction": net,
              "missing_species": [key for key in keys if key not in selected]}
    if status["missing_species"]:
        status["status"] = "missing_eligible_certificate"
        return [], status
    balances = {name: reaction_balance([selected[key] for key in lhs], [selected[key] for key in rhs])
                for name, (lhs, rhs) in {**steps, "net_reaction": net}.items()}
    cancellation = species_stoichiometry(steps.values()) == species_stoichiometry([net])
    status.update(balances=balances, species_cancellation_exact=cancellation)
    if not cancellation or not all(value["balanced"] for value in balances.values()):
        status["status"] = "unbalanced_cycle_rejected"
        return [], status
    sources = {key: {"certificate": relative(selected[key].certificate),
                     "certificate_sha256": selected[key].ledger["certificate_sha256"],
                     "geometry_sha256": selected[key].ledger["accepted_geometry_sha256"],
                     "spectrum_sha256": selected[key].ledger["spectrum_sha256"]} for key in keys}
    rows = []
    for temperature in DEFAULT_TEMPERATURES:
        deltas = {name: reaction_delta([selected[key] for key in lhs], [selected[key] for key in rhs], temperature)
                  for name, (lhs, rhs) in steps.items()}
        net_delta = reaction_delta([selected[key] for key in net[0]], [selected[key] for key in net[1]], temperature)
        total = math.fsum(deltas.values())
        residual = total - net_delta
        rows.append({"catalyst_id": identifier, "temperature_K": temperature,
                     **{f"{name}_delta_G_kcal_mol": value for name, value in deltas.items()},
                     "cycle_sum_delta_G_kcal_mol": total, "net_reaction_delta_G_kcal_mol": net_delta,
                     "closure_residual_kcal_mol": residual, "closure_tolerance_kcal_mol": CYCLE_TOLERANCE_KCAL_MOL,
                     "closure_passed": abs(residual) <= CYCLE_TOLERANCE_KCAL_MOL,
                     "atoms_balanced_every_step": True, "charge_balanced_every_step": True,
                     "species_cancellation_exact": cancellation,
                     "net_equation": "benzyl_alcohol + aniline -> product_amine + water",
                     "source_certificates": json.dumps(sources, sort_keys=True),
                     "scope": "State-function accounting from the same actual minimum certificates; closure is algebraic consistency, not a verified catalytic pathway, TS, rate or TOF"})
    status.update(status="exported" if all(row["closure_passed"] for row in rows) else "cycle_energy_closure_failed",
                  maximum_absolute_residual_kcal_mol=max(abs(row["closure_residual_kcal_mol"]) for row in rows))
    return rows, status


def write_csv(path, rows, empty_fields=("status",)):
    fields = list(rows[0]) if rows else list(empty_fields)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def export(output):
    output.mkdir(parents=True, exist_ok=True)
    history, snapshots = {}, []
    for name in ("data/campaign_alpb_toluene/records.jsonl", "data/targeted_rescue/records.jsonl"):
        path = REPO / name
        if path.exists():
            raw = path.read_bytes()
            snapshots.append({"path": name, "snapshot_sha256": hashlib.sha256(raw).hexdigest(), "snapshot_bytes": len(raw)})
            lines = raw.decode("utf-8").splitlines()
            for line_index, line in enumerate(lines):
                if line.strip():
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError:
                        # An append-only live history can end in an unfinished write.
                        # Never ignore corruption in a complete line.
                        if line_index != len(lines)-1 or raw.endswith((b"\n", b"\r")):
                            raise
                        snapshots[-1]["incomplete_final_line_excluded"] = True
                        continue
                    native = record.get("native_result", {})
                    if native.get("output_directory") and record.get("metadata"):
                        history[historical_key(native["output_directory"])] = record["metadata"]
    ledger, certificates = [], []
    jobs = [(p, False) for root in ROOTS for p in sorted((REPO / root).glob("*__*/result.json"))]
    jobs += [(p, True) for p in sorted((REPO / "data/reference_thermochemistry").glob("*/stationarity/result.json"))]
    for path, reference in jobs:
        entry, certificate = validate_certificate(path, history, reference=reference)
        ledger.append(entry)
        if certificate is not None:
            certificates.append(certificate)
    selected = select_certificates(certificates)
    thermal_rows, reference_rows, states, reactions, reaction_status = [], [], [], [], []
    cycle_rows, cycle_status = [], []
    for key, certificate in sorted(selected.items()):
        source = certificate.data["source_metadata"]
        for temperature, row in sorted(certificate.rows.items()):
            exported = {"species_key": key, "catalyst_id": certificate.metadata.get("catalyst_id", ""),
                "state": certificate.metadata.get("state", certificate.kind), "temperature_K": temperature,
                **{field: row[field] for field in ("E_elec", "ZPVE", "H_298", "S_harmonic", "S_qRRHO", "G_298_qRRHO_sol", "delta_G_concentration")},
                "energy_unit": "kcal/mol", "entropy_unit": "cal/(mol K)", "charge": certificate.charge,
                "formula": certificate.atoms.get_chemical_formula(), "method": "GFN2-xTB/ALPB toluene gsolv + qRRHO + 1M",
                "certificate": relative(certificate.certificate), "certificate_sha256": certificate.ledger["certificate_sha256"],
                "native_protocol_reference": f"{relative(certificate.certificate)}#/attempts/{certificate.data['accepted_attempt']}/native_command",
                "export_protocol_reference": relative(output / "export_protocol.json"),
                "geometry": relative(certificate.geometry), "geometry_sha256": certificate.ledger["accepted_geometry_sha256"],
                "spectrum": relative(certificate.spectrum), "spectrum_sha256": certificate.ledger["spectrum_sha256"],
                "conformer_source": source.get("conformer", PurePosixPath(source.get("optimization_output", "").replace("\\", "/")).parent.name),
                "accepted_attempt": certificate.data["accepted_attempt"], "selection_temperature_K": BASE_TEMPERATURE,
                "selection_scope": "lowest G383.15 among eligible certificates; same geometry at all T",
                "pressure_pa": row["pressure_pa"], "concentration_mol_l": row["concentration_mol_l"], "cutoff_cm1": row["cutoff_cm1"],
                "unpaired_electrons": certificate.data["unpaired_electrons"], "symmetry_number": certificate.data["symmetry_number"]}
            (reference_rows if certificate.kind == "reference" else thermal_rows).append(exported)
    schemes = ("alcohol_dehydrogenation", "separated_ion_deprotonation", "imine_hydrogenation")
    for spec in enumerate_specs():
        identifier = spec.catalyst_id
        for state in STATES:
            key = f"{identifier}|{state}"
            selected_certificate = selected.get(key)
            states.append({"catalyst_id": identifier, "state": state, "status": "certified" if selected_certificate else "no_eligible_certificate",
                "certificate": relative(selected_certificate.certificate) if selected_certificate else "",
                "reasons": json.dumps([entry.get("reason", entry["validation_status"]) for entry in ledger if entry["key"] == key and not entry["eligible"]])})
        for scheme in schemes:
            if scheme == "alcohol_dehydrogenation":
                lhs, rhs = [f"{identifier}|active", "benzyl_alcohol"], [f"{identifier}|hydrogenated", "benzaldehyde"]
                equation = "active_cat + benzyl_alcohol -> hydrogenated_cat + benzaldehyde"
                scope = "Reaction thermodynamics of separate 1M molecular species; no kinetic barrier or TOF"
            elif scheme == "separated_ion_deprotonation":
                lhs, rhs = [f"{identifier}|protonated_reference", "tert_butoxide_anion"], [f"{identifier}|active", "tert_butanol"]
                equation = "protonated_cat(+) + tert_butoxide(-) -> active_cat + tert_butanol"
                scope = "Separated-ion reference only; does not establish dissolved tBuOK salt speciation or ion pairing"
            else:
                lhs, rhs = [f"{identifier}|hydrogenated", "imine"], [f"{identifier}|active", "product_amine"]
                equation = "hydrogenated_cat + imine -> active_cat + product_amine"
                scope = "Reaction thermodynamics of separate 1M molecular species; no kinetic barrier or TOF"
            missing = [key for key in lhs + rhs if key not in selected]
            status = {"catalyst_id": identifier, "reaction": scheme, "equation": equation, "missing_species": missing}
            if missing:
                status["status"] = "missing_eligible_certificate"
            else:
                reactants, products = [selected[key] for key in lhs], [selected[key] for key in rhs]
                balance = reaction_balance(reactants, products)
                status.update(status="exported" if balance["balanced"] else "unbalanced_reaction_rejected", balance=balance)
                if balance["balanced"]:
                    for temperature in DEFAULT_TEMPERATURES:
                        delta = reaction_delta(reactants, products, temperature)
                        reactions.append({"catalyst_id": identifier, "reaction": scheme, "temperature_K": temperature,
                            "delta_G_kcal_mol": delta, "equation": equation, "atoms_balanced": True, "charge_balanced": True,
                            "reactant_certificates": json.dumps([relative(c.certificate) for c in reactants]),
                            "product_certificates": json.dumps([relative(c.certificate) for c in products]),
                            "reactant_geometry_sha256": json.dumps([c.ledger["accepted_geometry_sha256"] for c in reactants]),
                            "product_geometry_sha256": json.dumps([c.ledger["accepted_geometry_sha256"] for c in products]), "scope": scope})
            reaction_status.append(status)
        rows_for_cycle, status_for_cycle = export_cycle(selected, identifier)
        cycle_rows.extend(rows_for_cycle)
        cycle_status.append(status_for_cycle)
    write_csv(output / "catalyst_thermochemistry.csv", thermal_rows)
    write_csv(output / "reference_thermochemistry.csv", reference_rows)
    write_csv(output / "design_state_status.csv", states)
    write_csv(output / "reaction_free_energies.csv", reactions)
    write_csv(output / "thermodynamic_cycle_closure.csv", cycle_rows)
    atomic_json(output / "certificates.json", ledger)
    atomic_json(output / "selected_certificates.json", {key: certificate.ledger for key, certificate in sorted(selected.items())})
    atomic_json(output / "reaction_status.json", reaction_status)
    atomic_json(output / "cycle_status.json", cycle_status)
    summary = {"created_utc": utc_now(), "selection_temperature_K": BASE_TEMPERATURE, "temperature_grid_K": DEFAULT_TEMPERATURES,
        "selection": "Lowest certified G at 383.15 K per species. This is best among certified minima, not proof of a global minimum or the latest electronic-energy best. The selected geometry is fixed across all 13 temperatures.",
        "electronic_model": "Actual GFN2-xTB/ALPB toluene gsolv, nominal occupation and supplied symmetry, 300 K electronic smearing",
        "thermal_model": "Full projected molecular Hessian; ASE molecular RRHO with Grimme entropy-only qRRHO cutoff100cm-1; one ideal 1atm-to-1M correction",
        "limitations": "Solvent temperature derivatives, ion-pair equilibria, ground spin, global minima, kinetic barriers and TOF are not established.",
        "validation": "Local geometry/NPZ/native-output SHA256, actual gradient trace, recomputed full spectrum and thermochemistry, all covalent/ancillary identity and donor-coordination screens; raw records unchanged",
        "numerical_reproduction_tolerance": "1e-4 in source kcal/mol or cal/(mol K); Hessian frequency tolerance 1e-4cm-1 + 1e-7 relative; XYZ/NPZ position tolerance1e-7A",
        "history_snapshots": snapshots, "certificate_count": len(ledger), "eligible_certificate_count": len(certificates),
        "selected_catalyst_states": sum(c.kind == "catalyst" for c in selected.values()),
        "selected_active_catalysts": sum(c.kind == "catalyst" and c.metadata["state"] == "active" for c in selected.values()),
        "selected_reference_species": sum(c.kind == "reference" for c in selected.values()),
        "catalyst_temperature_rows": len(thermal_rows), "reference_temperature_rows": len(reference_rows), "reaction_temperature_rows": len(reactions),
        "exported_reactions": sum(row["status"] == "exported" for row in reaction_status),
        "reaction_schemes": schemes,
        "cycle_temperature_rows": len(cycle_rows),
        "exported_cycles": sum(row["status"] == "exported" for row in cycle_status),
        "cycle_closure_tolerance_kcal_mol": CYCLE_TOLERANCE_KCAL_MOL,
        "cycle_maximum_absolute_residual_kcal_mol": max((abs(row["closure_residual_kcal_mol"]) for row in cycle_rows), default=None),
        "cycle_interpretation": "Four component reaction free energies telescope to the net alcohol/aniline-to-amine/water reaction using identical species certificates. This is arithmetic/state-function consistency, not pathway or kinetic validation.",
        "files": {path.name: digest(path) for path in output.iterdir() if path.is_file() and path.name != "export_protocol.json"}}
    atomic_json(output / "export_protocol.json", summary)
    print(json.dumps({key: summary[key] for key in ("certificate_count", "eligible_certificate_count", "selected_catalyst_states", "selected_active_catalysts", "selected_reference_species", "reaction_temperature_rows")}))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPO / "data/datasets/thermochemistry")
    args = parser.parse_args()
    export(args.output.resolve())
