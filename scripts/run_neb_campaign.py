"""Real GFN2-xTB path searches from atom-mapped alcohol/catalyst endpoints.

Geometric endpoint construction is a hypothesis. Each complete endpoint must
subsequently optimize and retain the annotated chemical basin before CI-NEB.
This driver never reports a free-energy barrier or fabricates missing energies.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import gzip
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from ase import Atoms
from ase.constraints import FixBondLengths
from ase.data import covalent_radii
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import FIRE
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from pincer_catmech.kinetics.neb_ts_search import run_path
from pincer_catmech.quantum.xtb_backend import XTBCalculator, atomic_json, optimize_geometry

REPO = Path(__file__).resolve().parents[1]
DEFAULT_CATALYSTS = tuple(f"{metal}_macho_pnp_{substituent}"
                         for metal in ("Ru", "Mn", "Co") for substituent in ("Ph", "iPr"))


@dataclass
class ReactionAssembly:
    reactant: Atoms
    product: Atoms
    transfer_indices: dict[str, int]
    metadata: dict


def _utc():
    return datetime.now(timezone.utc).isoformat()


def ligand_identity_screen(atoms, metadata):
    """Audit the original ligand and ancillaries, allowing mapped substrate H uptake."""
    from pincer_catmech.quantum.identity import catalyst_identity_screen
    return catalyst_identity_screen(atoms, metadata)


def load_restart_hessian_mode(candidate, source, transfer_indices):
    """Select an actual negative mode by transfer overlap, without certifying a TS."""
    from pincer_catmech.kinetics.neb_ts_search import _project_cartesian_hessian, transfer_mode_overlap
    source = Path(source).resolve()
    with np.load(source, allow_pickle=False) as archive:
        if (not np.array_equal(archive["atomic_numbers"], candidate.numbers)
                or not np.array_equal(archive["masses"], candidate.get_masses())
                or archive["positions"].shape != candidate.positions.shape
                or not np.allclose(archive["positions"], candidate.positions, rtol=0, atol=1e-9)):
            raise ValueError("Hessian source does not match the exact candidate atoms, masses and coordinates")
        matrix = archive["hessian_ev_a2"].copy()
        step = float(archive["displacement_angstrom"])
    if matrix.shape != (3 * len(candidate),) * 2 or not np.isfinite(matrix).all():
        raise ValueError("Restart mode requires a finite full Cartesian Hessian")
    projected = _project_cartesian_hessian(candidate, matrix, step)
    if projected.antisymmetry_relative > .05:
        raise ValueError("Restart Hessian antisymmetry exceeds the existing Hessian tolerance")
    rows = []
    for index in np.flatnonzero(projected.frequencies_cm1 < 0):
        overlaps = transfer_mode_overlap(candidate, projected.cartesian_modes[index], transfer_indices)
        rows.append({"index": int(index), "frequency_cm1": float(projected.frequencies_cm1[index]), **overlaps})
    if not rows:
        raise ValueError("Hessian contains no actual negative internal mode for saddle initialization")
    selected = max(rows, key=lambda row: row["combined_overlap"])
    return projected.cartesian_modes[selected["index"]].copy(), {
        "source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "selection": "Largest combined mapped proton/hydride overlap among independently reprojected negative internal modes",
        "negative_modes": rows, "selected_mode": selected,
        "scope": "Dimer initialization only; source may be nonstationary and fail every final TS gate"}


def make_native_hessian_provider(scratch_root, label, *, charge, unpaired,
                                threads, executable, deadline, solvent, solvation_state):
    """Use one native process for the full unbiased Hessian, retaining its evidence."""
    def provider(candidate, attempt_directory):
        from pincer_catmech.quantum.thermochemistry import native_hessian
        folder = Path(scratch_root) / label / "native_hessian" / Path(attempt_directory).name
        def archive():
            artifacts = {}
            for name in ("input.xyz", "hessian.inp", "command.json", "xtb.out", "hessian.out"):
                source = folder / name
                if source.exists():
                    data = source.read_bytes()
                    filename = "native_" + name + (".gz" if name in ("xtb.out", "hessian.out") else "")
                    saved = gzip.compress(data, mtime=0) if filename.endswith(".gz") else data
                    (Path(attempt_directory) / filename).write_bytes(saved)
                    artifacts[name] = {"file": filename, "source_sha256": hashlib.sha256(data).hexdigest(),
                                       "archive_sha256": hashlib.sha256(saved).hexdigest()}
            return artifacts
        gradient_energy = float(candidate.get_potential_energy())
        try:
            result, native_energy, elapsed = native_hessian(candidate, folder,
                charge=charge, unpaired=unpaired, threads=threads, executable=executable,
                deadline_epoch=deadline, timeout_seconds=900,
                solvent=solvent, solvation_state=solvation_state)
        except Exception as error:
            atomic_json(Path(attempt_directory) / "native_hessian_provider_failure.json",
                        {"error": f"{type(error).__name__}: {error}", "native_archives": archive()})
            raise
        difference = float(native_energy - gradient_energy)
        evidence = {"provider": "native unbiased GFN2-xTB numerical full Hessian",
            "scratch_directory": str(folder), "native_electronic_energy_eV": native_energy,
            "gradient_electronic_energy_eV": gradient_energy, "energy_difference_eV": difference,
            "energy_consistency_tolerance_eV": 0.001, "wall_seconds": elapsed,
            "full_matrix_shape": list(result.hessian_ev_a2.shape),
            "finite_difference_step_A": result.displacement_angstrom,
            "solvent": solvent, "solvation_state": solvation_state,
            "required_native_checks": "Normal termination; full unbiased Hessian; zero frozen atoms; scale 1.0; finite full (3N)^2 unprojected matrix",
            "native_archives": archive(),
            "native_sha256": {name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                              for name in ("input.xyz", "hessian.inp", "command.json", "xtb.out", "hessian.out")}}
        atomic_json(Path(attempt_directory) / "native_hessian_provider.json", evidence)
        if abs(difference) > 0.001:
            raise ValueError("Native Hessian electronic energy differs from the gradient surface by more than 0.001 eV")
        return result
    return provider


def _benzyl_alcohol(seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("[OH:1][CH2:2]c1ccccc1"))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(seed) % (2**31 - 1)
    if AllChem.EmbedMolecule(mol, parameters) != 0:
        raise RuntimeError("RDKit failed to embed benzyl alcohol")
    AllChem.UFFOptimizeMolecule(mol, maxIters=300)
    mapped = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    oxygen, carbon = mapped[1], mapped[2]
    proton = next(a.GetIdx() for a in mol.GetAtomWithIdx(oxygen).GetNeighbors() if a.GetAtomicNum() == 1)
    hydrides = [a.GetIdx() for a in mol.GetAtomWithIdx(carbon).GetNeighbors() if a.GetAtomicNum() == 1]
    return mol, oxygen, carbon, proton, hydrides


def _place_substrate(catalyst, mol, oxygen, carbon, proton, hydrides, site, seed):
    """Collision-aware rigid placement using geometry, with no energy surrogate."""
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    local = mol.GetConformer().GetPositions() - mol.GetConformer().GetPositions()[carbon]
    metal_position, site_position = catalyst.positions[0], catalyst.positions[site]
    thresholds = 1.03 * (covalent_radii[catalyst.numbers, None] + covalent_radii[numbers][None, :])
    thresholds = np.maximum(thresholds, .85)
    rng = np.random.default_rng(seed)
    best = None
    for hydride in hydrides:
        def coordinates(parameters):
            return local @ Rotation.from_rotvec(parameters[:3]).as_matrix().T + parameters[3:]

        def score(parameters):
            placed = coordinates(parameters)
            distances = np.linalg.norm(catalyst.positions[:, None, :] - placed[None, :, :], axis=2)
            overlap = np.maximum(thresholds - distances, 0)
            direction_h = placed[carbon] - metal_position
            direction_p = placed[oxygen] - site_position
            endpoint_h = metal_position + 1.60 * direction_h / max(np.linalg.norm(direction_h), 1e-8)
            endpoint_p = site_position + 1.02 * direction_p / max(np.linalg.norm(direction_p), 1e-8)
            projected = np.stack((endpoint_h, endpoint_p))
            projected_distances = np.linalg.norm(catalyst.positions[:, None, :] - projected[None, :, :], axis=2)
            projected_thresholds = np.maximum(covalent_radii[catalyst.numbers, None] + .31, 1.0)
            projected_overlap = np.maximum(projected_thresholds - projected_distances, 0)
            projected_overlap[0, 0] = 0
            projected_overlap[site, 1] = 0
            return float(200 * np.sum(overlap**2)
                + 300 * np.sum(projected_overlap**2)
                + 8 * (np.linalg.norm(placed[hydride] - metal_position) - 1.80)**2
                + 8 * (np.linalg.norm(placed[proton] - site_position) - 1.75)**2
                + .5 * (np.linalg.norm(placed[oxygen] - site_position) - 2.70)**2)

        candidates = []
        for _ in range(160):
            rotation = Rotation.random(random_state=rng).as_rotvec()
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            translation = metal_position + 2.85 * direction
            parameters = np.concatenate((rotation, translation))
            candidates.append((score(parameters), parameters))
        for _, initial in sorted(candidates, key=lambda row: row[0])[:5]:
            optimized = minimize(score, initial, method="L-BFGS-B", options={"maxiter": 220, "ftol": 1e-10})
            placed = coordinates(optimized.x)
            minimum = float(np.linalg.norm(catalyst.positions[:, None, :] - placed[None, :, :], axis=2).min())
            if minimum > .80 and (best is None or optimized.fun < best[0]):
                best = (float(optimized.fun), placed.copy(), hydride, minimum)
    if best is None:
        raise ValueError("No collision-free alcohol encounter pose was found")
    return best


def _aldehyde_coordinates(mol, coordinates, oxygen, carbon, proton, hydride, seed):
    """Embed the exact mapped neutral aldehyde graph and align its retained atoms."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    editable.GetBondBetweenAtoms(oxygen, carbon).SetBondType(Chem.BondType.DOUBLE)
    for index in sorted((proton, hydride), reverse=True):
        editable.RemoveAtom(index)
    aldehyde = editable.GetMol()
    Chem.SanitizeMol(aldehyde)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = (int(seed) + 117) % (2**31 - 1)
    if AllChem.EmbedMolecule(aldehyde, parameters) != 0:
        raise RuntimeError("RDKit failed to embed the mapped benzaldehyde product")
    AllChem.UFFOptimizeMolecule(aldehyde, maxIters=300)
    original_indices = np.array([a.GetAtomMapNum() - 1 for a in aldehyde.GetAtoms()])
    source = aldehyde.GetConformer().GetPositions()
    heavy = np.array([a.GetAtomicNum() > 1 for a in aldehyde.GetAtoms()])
    target = coordinates[original_indices]
    center_a, center_b = source[heavy].mean(0), target[heavy].mean(0)
    u, _, vt = np.linalg.svd((source[heavy] - center_a).T @ (target[heavy] - center_b))
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    aligned = (source - center_a) @ (u @ vt) + center_b
    return original_indices, aligned


def assemble_endpoints(catalyst: Atoms, structure_metadata: dict, *, seed: int = 4097) -> ReactionAssembly:
    """Construct full Cat+PhCH2OH and Cat-H2+PhCHO endpoints with identical atoms."""
    if np.any(catalyst.pbc) or catalyst.constraints or not np.isfinite(catalyst.positions).all():
        raise ValueError("Catalyst must be a finite unconstrained molecule")
    if structure_metadata.get("state") != "active":
        raise ValueError("The reaction requires the explicitly deprotonated active state")
    site, metal = int(structure_metadata["proton_site_index"]), int(structure_metadata["metal_index"])
    if metal != 0 or not 0 <= site < len(catalyst):
        raise ValueError("Invalid mapped metal/proton site")
    if catalyst[site].symbol != structure_metadata["proton_site_element"]:
        raise ValueError("Catalyst atom ordering does not match proton-site metadata")
    mol, oxygen, carbon, proton, hydrides = _benzyl_alcohol(seed)
    score, placed, hydride, minimum = _place_substrate(catalyst, mol, oxygen, carbon, proton, hydrides, site, seed)
    substrate = Atoms([a.GetSymbol() for a in mol.GetAtoms()], positions=placed)
    reactant = catalyst.copy() + substrate
    product = reactant.copy()
    offset = len(catalyst)
    retained, aldehyde = _aldehyde_coordinates(mol, placed, oxygen, carbon, proton, hydride, seed)
    product.positions[offset + retained] = aldehyde
    site_position, metal_position = product.positions[site], product.positions[metal]
    proton_direction = placed[oxygen] - site_position
    proton_direction /= np.linalg.norm(proton_direction)
    proton_length = {"N": 1.02, "C": 1.09, "O": .98}[catalyst[site].symbol]
    product.positions[offset + proton] = site_position + proton_length * proton_direction
    hydride_direction = placed[carbon] - metal_position
    hydride_direction /= np.linalg.norm(hydride_direction)
    product.positions[offset + hydride] = metal_position + 1.60 * hydride_direction
    # Move the retained aldehyde as a rigid fragment until both transferred
    # hydrogens belong geometrically to the catalyst, avoiding a disguised
    # reactant basin. Subsequent unconstrained quantum optimization decides
    # whether this encounter arrangement is an actual local minimum.
    departure = proton_direction + hydride_direction
    departure_length = np.linalg.norm(departure)
    if departure_length < 1e-8:
        raise ValueError("Opposed transfer directions do not define a product departure direction")
    departure /= departure_length
    departure_distance = None
    for displacement in (0.0, .35, .70, 1.05, 1.40, 1.80, 2.20, 2.80, 3.40):
        product.positions[offset + retained] = aldehyde + displacement * departure
        proton_progress = (product.get_distance(offset + proton, offset + oxygen)
                           - product.get_distance(offset + proton, site))
        hydride_progress = (product.get_distance(offset + hydride, offset + carbon)
                            - product.get_distance(offset + hydride, metal))
        minimum_product = np.min(product.get_all_distances() + np.eye(len(product)) * 999)
        if proton_progress > .30 and hydride_progress > .30 and minimum_product > .75:
            departure_distance = displacement
            break
    if departure_distance is None:
        raise ValueError("Constructed product cannot bracket both mapped hydrogen transfers")
    electronic = {"charge": int(structure_metadata["charge"]),
                  "unpaired": int(structure_metadata["unpaired_electrons"]),
                  "uhf": int(structure_metadata["unpaired_electrons"])}
    reactant.info = {**electronic, "endpoint_role": "alcohol encounter complex; requires optimization"}
    product.info = {**electronic, "endpoint_role": "hydrogenated catalyst plus aldehyde; requires optimization"}
    transfer = {"proton": offset + proton, "proton_donor": offset + oxygen, "proton_acceptor": site,
                "hydride": offset + hydride, "hydride_donor": offset + carbon, "hydride_acceptor": metal}
    distances = {}
    for endpoint_name, endpoint in (("reactant", reactant), ("product", product)):
        pair_distances = endpoint.get_all_distances() + np.eye(len(endpoint)) * 999
        if pair_distances.min() < .70:
            i, j = np.unravel_index(np.argmin(pair_distances), pair_distances.shape)
            raise ValueError(f"Constructed {endpoint_name} has an atom overlap below 0.70 angstrom: "
                             f"{i}({endpoint[i].symbol})-{j}({endpoint[j].symbol})={pair_distances[i,j]:.4f}")
        distances[endpoint_name] = {
            channel: {"donor_H_A": float(endpoint.get_distance(transfer[channel + "_donor"], transfer[channel])),
                      "acceptor_H_A": float(endpoint.get_distance(transfer[channel + "_acceptor"], transfer[channel]))}
            for channel in ("proton", "hydride")}
    metadata = {"seed": seed, "catalyst_id": structure_metadata["catalyst_id"],
        "catalyst_atom_count": offset, "complete_atom_count": len(reactant),
        "reactant_formula": reactant.get_chemical_formula(), "product_formula": product.get_chemical_formula(),
        "transfer_indices": transfer, "substrate_oxygen_index": offset + oxygen,
        "substrate_carbonyl_carbon_index": offset + carbon, "electronic_state": electronic,
        "catalyst_metadata": structure_metadata,
        "placement_score": score, "placement_score_units": "geometric penalty; NOT an energy",
        "minimum_initial_interfragment_distance_A": minimum, "input_transfer_distances": distances,
        "aldehyde_rigid_departure_A": departure_distance,
        "construction_status": "Mapped geometric hypotheses; endpoint optimization and basin validation are mandatory",
        "substrate_preparation": "RDKit ETKDG plus organic UFF geometry preparation; no UFF energies used as quantum data"}
    return ReactionAssembly(reactant, product, transfer, metadata)


def _load_best(path):
    if not Path(path).exists():
        return {}
    records = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    return {record["catalyst_id"]: record for record in records.values()
            if record.get("state") == "active" and record["native_result"]["converged"] and record.get("coordination_retained")}


def run_one(record, *, output_root, scratch_root, executable, threads, deadline_epoch,
            search_seconds, pilot=False, seed=4097, solvent="toluene", solvation_state="gsolv",
            restrained_product_preconditioning=False, refine_unconverged_band=False,
            hessian_method="ase"):
    identifier = record["catalyst_id"]
    label = f"{identifier}_{'pilot' if pilot else 'search'}_{int(time.time())}_{os.getpid()}"
    directory = Path(output_root) / label
    directory.mkdir(parents=True, exist_ok=False)
    source = Path(record["native_result"]["output_directory"]) / "xtbopt.xyz"
    started = time.time()
    task = {"catalyst_id": identifier, "mode": "pilot" if pilot else "search", "started_utc": _utc(),
            "pid": os.getpid(), "output_directory": str(directory), "source_geometry": str(source),
            "source_geometry_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "source_conformer": record["conformer"], "threads": threads,
            "solvent": solvent, "solvation_model": "ALPB", "solvation_state": solvation_state,
            "restrained_product_preconditioning_requested": restrained_product_preconditioning,
            "refine_unconverged_band": refine_unconverged_band,
            "hessian_method": hessian_method,
            "source_geometry_role": "Preoptimization only; both complete endpoints reevaluated under the declared solvent convention",
            "no_free_energy_barrier": True, "status": "assembling"}
    atomic_json(directory / "task.json", task)
    try:
        expected_hash = record.get("optimized_geometry_sha256")
        if expected_hash is not None and expected_hash != task["source_geometry_sha256"]:
            raise ValueError("Source optimized geometry hash differs from the accepted conformer record")
        task["source_ligand_identity"] = ligand_identity_screen(read(source), record["metadata"])
        if not task["source_ligand_identity"]["retained"]:
            raise ValueError("Source active structure changed its ligand covalent identity; select an identity-preserving conformer")
        assembly_errors = []
        assembly = None
        for attempt in range(4):
            try:
                assembly = assemble_endpoints(read(source), record["metadata"], seed=seed + attempt * 1031)
                break
            except ValueError as error:
                assembly_errors.append({"seed": seed + attempt * 1031, "error": str(error)})
        task["assembly_retries"] = assembly_errors
        if assembly is None:
            raise ValueError("All bounded collision-aware endpoint assembly attempts failed")
        write(directory / "constructed_reactant.xyz", assembly.reactant, format="xyz")
        write(directory / "constructed_product.xyz", assembly.product, format="xyz")
        atomic_json(directory / "assembly.json", assembly.metadata)
        allowed = min(search_seconds, deadline_epoch - time.time())
        if allowed <= 0:
            raise TimeoutError("Campaign deadline reached before the path search")
        search_deadline = min(deadline_epoch, time.time() + allowed)
        if not pilot:
            task["endpoint_preoptimization"] = []
            for index, endpoint in enumerate((assembly.reactant, assembly.product)):
                remaining = search_deadline - time.time()
                if remaining <= 0:
                    raise TimeoutError("Search deadline reached during endpoint preoptimization")
                if index == 1 and restrained_product_preconditioning:
                    transfer = assembly.transfer_indices
                    pairs = [(transfer["hydride_acceptor"], transfer["hydride"]),
                             (transfer["proton_acceptor"], transfer["proton"])]
                    endpoint.calc = XTBCalculator(Path(scratch_root) / label / "product_preconditioning",
                        charge=record["metadata"]["charge"], unpaired=record["metadata"]["unpaired_electrons"],
                        threads=threads, executable=executable, timeout_seconds=120,
                        deadline=min(search_deadline, time.time() + 300), solvent=solvent,
                        solvation_state=solvation_state,
                        trace_path=directory / "evaluations_product_preconditioning.jsonl.gz")
                    constraint_record = {"method": "ASE FixBondLengths plus FIRE on real ALPB GFN2-xTB potential",
                        "fixed_pairs_zero_based": pairs, "bond_lengths_A": [float(endpoint.get_distance(*pair)) for pair in pairs],
                        "energy_use": "Preparation only; excluded from accepted barrier/free-energy results"}
                    endpoint.set_constraint(FixBondLengths(pairs))
                    try:
                        with FIRE(endpoint, logfile=str(directory / "product_preconditioning.log"),
                                  trajectory=str(directory / "product_preconditioning.traj")) as optimizer:
                            constraint_record["projected_force_converged"] = bool(optimizer.run(fmax=.08, steps=180))
                            constraint_record["steps"] = optimizer.nsteps
                    finally:
                        endpoint.set_constraint()
                        constraint_record["constraints_removed_before_validation"] = not endpoint.constraints
                        atomic_json(directory / "product_preconditioning.json", constraint_record)
                    write(directory / "preconditioned_product.xyz", endpoint, format="xyz")
                    endpoint.calc = None
                remaining = search_deadline - time.time()
                if remaining <= 0:
                    raise TimeoutError("Search deadline reached before unconstrained endpoint optimization")
                native_directory = directory / f"native_endpoint_{index}"
                native = optimize_geometry(endpoint, native_directory,
                    charge=record["metadata"]["charge"], unpaired=record["metadata"]["unpaired_electrons"],
                    threads=threads, executable=executable, max_cycles=250,
                    timeout_seconds=min(300, remaining), solvent=solvent, solvation_state=solvation_state)
                task["endpoint_preoptimization"].append(asdict(native))
                print(json.dumps({"utc": _utc(), "catalyst_id": identifier,
                    "stage": f"native_endpoint_{index}", "converged": native.converged,
                    "energy_eV": native.energy_eV, "max_force_eV_A": native.max_force_eV_A}), flush=True)
                saved_native_geometry = native_directory / "xtbopt.xyz"
                if saved_native_geometry.exists():
                    optimized = read(saved_native_geometry)
                    if not np.array_equal(optimized.numbers, endpoint.numbers):
                        raise ValueError("Native endpoint optimization changed atomic identity/order")
                    endpoint_record = task["endpoint_preoptimization"][-1]
                    endpoint_record["ligand_identity"] = ligand_identity_screen(optimized, record["metadata"])
                    if native.converged and not endpoint_record["ligand_identity"]["retained"]:
                        raise ValueError("Optimized endpoint changed the original ligand covalent identity")
                    if endpoint_record["ligand_identity"]["retained"]:
                        endpoint.positions = optimized.positions
                        endpoint_record["used_geometry_sha256"] = hashlib.sha256(saved_native_geometry.read_bytes()).hexdigest()
                        endpoint_record["used_geometry_role"] = ("Converged native seed; fresh force/basin checks remain mandatory"
                            if native.converged else "Actual last unconverged native seed; no convergence claim; fresh force/basin checks remain mandatory")
                    else:
                        endpoint_record["used_geometry_role"] = "Unconverged native seed changed chemical identity and was not used"
                atomic_json(directory / "task.json", task)
            allowed = search_deadline - time.time()
            if allowed <= 0:
                raise TimeoutError("Search deadline reached after endpoint preoptimization")

        def factory(calculator_label):
            return XTBCalculator(Path(scratch_root) / label / calculator_label,
                charge=record["metadata"]["charge"], unpaired=record["metadata"]["unpaired_electrons"],
                threads=threads, executable=executable, timeout_seconds=120,
                deadline=search_deadline, solvent=solvent, solvation_state=solvation_state,
                trace_path=directory / f"evaluations_{calculator_label}.jsonl.gz")

        last_emit = {"stage": None, "time": 0.0}

        def telemetry(event):
            if event["stage"] != last_emit["stage"] or time.time() - last_emit["time"] >= 600:
                print(json.dumps({"utc": _utc(), "catalyst_id": identifier, **event}, allow_nan=False), flush=True)
                last_emit.update(stage=event["stage"], time=time.time())

        result = run_path(assembly.reactant, assembly.product, factory, directory / "path",
            assembly.transfer_indices, endpoint_steps=8 if pilot else 200,
            neb_steps=4 if pilot else 120, ts_steps=4 if pilot else 80,
            max_corrections=0 if pilot else 2, descent_steps=80,
            refine_unconverged_band=refine_unconverged_band,
            hessian_provider=make_native_hessian_provider(scratch_root, label,
                charge=record["metadata"]["charge"], unpaired=record["metadata"]["unpaired_electrons"],
                threads=threads, executable=executable, deadline=search_deadline,
                solvent=solvent, solvation_state=solvation_state) if hessian_method == "native" else None,
            deadline_seconds=allowed, progress_callback=telemetry)
        task.update(status=result.status, accepted=result.accepted, path_result=asdict(result))
    except Exception as error:
        task.update(status="driver_failed", accepted=False, error=f"{type(error).__name__}: {error}")
    task.update(finished_utc=_utc(), elapsed_seconds=time.time() - started)
    atomic_json(directory / "task.json", task)
    return task


def resume_one(source_route, *, output_root, scratch_root, executable, threads,
               deadline_epoch, search_seconds=2700, solvent="toluene", solvation_state="gsolv",
               endpoints_only=False, ts_initial_image=None, resume_last_dimer=False,
               resume_neb_steps=None, hessian_method="ase", native_endpoints=False,
               restart_hessian=None):
    """Continue a saved band, or explicitly retry its two optimized endpoints."""
    source_route = Path(source_route).resolve()
    source_task = json.loads((source_route / "task.json").read_text(encoding="utf-8"))
    source_result = json.loads((source_route / "path/result.json").read_text(encoding="utf-8"))
    if source_result["accepted"]:
        raise ValueError("The source path is already accepted; no continuation is necessary")
    if source_task.get("solvent") != solvent or source_task.get("solvation_state") != solvation_state:
        raise ValueError("Continuation solvent convention must equal the source path convention")
    assembly = json.loads((source_route / "assembly.json").read_text(encoding="utf-8"))
    trajectory_path = source_route / "path/neb.traj"
    if endpoints_only:
        endpoint_files = [(source_route / f"native_endpoint_{index}/xtbopt.xyz" if native_endpoints
                           else source_route / f"path/endpoint_{index}.xyz") for index in (0, 1)]
        images = [read(path) for path in endpoint_files]
        restart_provenance = {"source_endpoint_sha256":
            [hashlib.sha256(path.read_bytes()).hexdigest() for path in endpoint_files],
            "source_endpoint_files": [str(path) for path in endpoint_files],
            "source_endpoint_role": "Actual final native coordinates, possibly unconverged; mandatory fresh force/optimization/basin gates"
                                    if native_endpoints else "Previously optimized path endpoints; fresh force and basin gates"}
    else:
        with Trajectory(trajectory_path, "r") as trajectory:
            total_frames = len(trajectory)
            complete_blocks = total_frames // 9
            if complete_blocks < 1:
                raise ValueError("Source trajectory has no complete nine-image band; explicitly select --resume-endpoints to retry saved endpoints")
            start = (complete_blocks - 1) * 9
            images = [trajectory[index] for index in range(start, start + 9)]
        restart_provenance = {
            "source_trajectory_sha256": hashlib.sha256(trajectory_path.read_bytes()).hexdigest(),
            "source_trajectory_frame_count": total_frames,
            "selected_frame_indices": list(range(start, start + 9))}
    for image in images:
        image.info = dict(assembly["electronic_state"])
        image.calc = None
    candidate = None
    candidate_source = None
    if resume_last_dimer:
        if endpoints_only:
            raise ValueError("A dimer candidate continuation requires its saved nine-image band")
        candidate_sources = sorted((source_route / "path").glob("attempt_*/dimer.traj"), reverse=True)
        for candidate_source in candidate_sources:
            with Trajectory(candidate_source, "r") as trajectory:
                if len(trajectory):
                    candidate = trajectory[-1]
                    candidate_frame = len(trajectory) - 1
                    break
        if candidate is None:
            raise ValueError("No complete saved dimer geometry is available for continuation")
        candidate.info = dict(assembly["electronic_state"])
        candidate.calc = None
        restart_provenance["dimer_restart"] = {
            "source": str(candidate_source), "source_frame": candidate_frame,
            "source_sha256": hashlib.sha256(candidate_source.read_bytes()).hexdigest(),
            "scope": "Retain actual candidate coordinates; optimizer history and dimer mode are reinitialized"}
        if ts_initial_image is None:
            ts_initial_image = source_result.get("selected_initial_image") or source_task.get("requested_ts_initial_image")
    restart_mode = None
    if restart_hessian is not None:
        if candidate is None:
            raise ValueError("A Hessian mode restart requires the exact saved dimer candidate")
        restart_mode, mode_provenance = load_restart_hessian_mode(candidate, restart_hessian, assembly["transfer_indices"])
        restart_provenance["hessian_mode_restart"] = mode_provenance
    identifier = source_task["catalyst_id"]
    mode = "retry_endpoints" if endpoints_only else "resume"
    label = f"{identifier}_{mode}_{int(time.time())}_{os.getpid()}"
    directory = Path(output_root) / label
    directory.mkdir(parents=True, exist_ok=False)
    remaining = min(search_seconds, deadline_epoch - time.time())
    if remaining <= 0:
        raise TimeoutError("Global campaign deadline reached before continuation")
    search_deadline = min(deadline_epoch, time.time() + remaining)
    task = {"catalyst_id": identifier, "mode": mode, "started_utc": _utc(), "pid": os.getpid(),
        "source_route": str(source_route), "source_result_status": source_result["status"],
        **restart_provenance,
        "source_protocol_sha256": hashlib.sha256((source_route / "path/protocol.json").read_bytes()).hexdigest(),
        "output_directory": str(directory), "threads": threads, "solvent": solvent,
        "solvation_model": "ALPB", "solvation_state": solvation_state, "no_free_energy_barrier": True,
        "refine_unconverged_band": True,
        "requested_ts_initial_image": ts_initial_image,
        "requested_resume_neb_steps": resume_neb_steps,
        "resume_last_dimer": resume_last_dimer,
        "hessian_method": hessian_method,
        "continuation_protocol": ("Reuse actual optimized endpoint coordinates; fresh force and transfer-basin gates; new IDPP band"
                                  if endpoints_only else
                                  "Reuse only actual saved coordinates; fresh forces on every band image; no native endpoint or IDPP repetition"),
        "status": "resuming"}
    atomic_json(directory / "task.json", task)
    atomic_json(directory / "assembly.json", assembly)
    write(directory / ("restart_endpoints.extxyz" if endpoints_only else "restart_band.extxyz"), images)

    def factory(calculator_label):
        return XTBCalculator(Path(scratch_root) / label / calculator_label,
            charge=assembly["electronic_state"]["charge"], unpaired=assembly["electronic_state"]["unpaired"],
            threads=threads, executable=executable, timeout_seconds=120, deadline=search_deadline,
            solvent=solvent, solvation_state=solvation_state,
            trace_path=directory / f"evaluations_{calculator_label}.jsonl.gz")

    try:
        catalyst_metadata = assembly.get("catalyst_metadata")
        if catalyst_metadata is None:
            from pincer_catmech.generators import CatalystSpec, generate_catalyst
            metal, remaining_name = identifier.split("_", 1)
            backbone, substituent = remaining_name.rsplit("_", 1)
            catalyst_metadata = generate_catalyst(CatalystSpec(metal, backbone, substituent), state="active").metadata()
        task["endpoint_ligand_identity"] = [ligand_identity_screen(endpoint, catalyst_metadata)
                                            for endpoint in (images[0], images[-1])]
        if not all(screen["retained"] for screen in task["endpoint_ligand_identity"]):
            raise ValueError("Saved endpoint changed the original ligand covalent identity; its path cannot be continued as the requested catalyst")
        result = run_path(images[0], images[-1], factory, directory / "path", assembly["transfer_indices"],
            restart_images=None if endpoints_only else images,
            restart_source=None if endpoints_only else trajectory_path,
            restart_ts_candidate=candidate, restart_ts_source=candidate_source,
            restart_ts_mode=restart_mode, restart_ts_mode_source=restart_hessian,
            endpoint_steps=80, neb_steps=resume_neb_steps if resume_neb_steps is not None else (120 if endpoints_only else 80), ts_steps=120,
            refine_unconverged_band=True, ts_initial_image=ts_initial_image,
            hessian_provider=make_native_hessian_provider(scratch_root, label,
                charge=assembly["electronic_state"]["charge"], unpaired=assembly["electronic_state"]["unpaired"],
                threads=threads, executable=executable, deadline=search_deadline,
                solvent=solvent, solvation_state=solvation_state) if hessian_method == "native" else None,
            max_corrections=2, descent_steps=100, deadline_seconds=remaining,
            progress_callback=lambda event: None)
        task.update(status=result.status, accepted=result.accepted, path_result=asdict(result))
    except Exception as error:
        task.update(status="driver_failed", accepted=False, error=f"{type(error).__name__}: {error}")
    task.update(finished_utc=_utc())
    atomic_json(directory / "task.json", task)
    return task


def run_resume_campaign(args, deadline_epoch):
    """Run at most two genuine checkpoint continuations without adding hidden workers."""
    scheduler = args.output_root / f"scheduler_resume_{int(time.time())}_{os.getpid()}.json"
    status = {"mode": "resume", "pid": os.getpid(), "started_utc": _utc(), "status": "running",
        "source_routes": [str(path.resolve()) for path in args.resume_routes], "workers": args.workers,
        "threads_per_worker": args.threads, "deadline_epoch": deadline_epoch, "results": []}
    atomic_json(scheduler, status)
    print(json.dumps(status), flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(resume_one, path, output_root=str(args.output_root), scratch_root=str(args.scratch_root),
            executable=args.xtb, threads=args.threads, deadline_epoch=deadline_epoch, search_seconds=args.search_seconds,
            solvent=args.solvent, solvation_state=args.solvation_state,
            endpoints_only=args.resume_endpoints or args.resume_native_endpoints,
            native_endpoints=args.resume_native_endpoints, ts_initial_image=args.ts_initial_image,
            resume_last_dimer=args.resume_last_dimer, resume_neb_steps=args.resume_neb_steps,
            hessian_method=args.hessian_method, restart_hessian=args.restart_hessian): str(path) for path in args.resume_routes}
        last_telemetry = time.time()
        while pending:
            completed, _ = wait(pending, timeout=2, return_when=FIRST_COMPLETED)
            for future in completed:
                source = pending.pop(future)
                try:
                    result = future.result()
                except Exception as error:
                    result = {"source_route": source, "status": "worker_failed", "accepted": False, "error": repr(error)}
                status["results"].append(result)
                print(json.dumps({"utc": _utc(), "completed_source": source, "status": result["status"]}), flush=True)
            status.update(active_source_routes=list(pending.values()), updated_utc=_utc())
            atomic_json(scheduler, status)
            if time.time() - last_telemetry >= 600:
                print(json.dumps({"utc": _utc(), "mode": "resume", "active_source_routes": list(pending.values())}), flush=True)
                last_telemetry = time.time()
    status.update(status="finished", finished_utc=_utc())
    atomic_json(scheduler, status)
    print(json.dumps(status), flush=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best", type=Path, default=REPO / "data/datasets/best_conformers.json")
    parser.add_argument("--control", type=Path, default=REPO / "data/campaign/control.json")
    parser.add_argument("--output-root", type=Path, default=REPO / "data/campaign/neb")
    parser.add_argument("--scratch-root", type=Path, default=REPO.parent / "neb-scratch")
    parser.add_argument("--xtb", default=os.environ.get("PINCER_XTB"))
    parser.add_argument("--catalysts", nargs="+", default=list(DEFAULT_CATALYSTS))
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--threads", type=int, choices=(1, 2), default=2)
    parser.add_argument("--search-seconds", type=float, default=1200)
    parser.add_argument("--wait-seconds", type=float, default=1800)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--solvent", choices=("toluene",), default="toluene")
    parser.add_argument("--solvation-state", choices=("gsolv",), default="gsolv")
    parser.add_argument("--seed", type=int, default=4097)
    parser.add_argument("--restrained-product-preconditioning", action="store_true")
    parser.add_argument("--resume-routes", type=Path, nargs="+")
    parser.add_argument("--resume-endpoints", action="store_true",
                        help="With --resume-routes, reuse only both actual optimized endpoints and initialize a fresh band")
    parser.add_argument("--resume-native-endpoints", action="store_true",
                        help="Reuse both actual native final coordinates, even if unconverged; fresh force/basin gates are mandatory")
    parser.add_argument("--ts-initial-image", type=int, choices=range(1, 8),
                        help="With --resume-routes, explicitly refine this local image; record the global peak separately")
    parser.add_argument("--resume-last-dimer", action="store_true",
                        help="Continue the last actual dimer geometry as well as its saved band; all TS gates remain mandatory")
    parser.add_argument("--resume-neb-steps", type=int,
                        help="Explicit additional NEB steps, e.g. 1 before resuming an actual dimer candidate")
    parser.add_argument("--restart-hessian", type=Path,
                        help="Use the exact candidate full Hessian's negative mode with largest mapped dual-transfer overlap")
    parser.add_argument("--refine-unconverged-band", action="store_true")
    parser.add_argument("--hessian-method", choices=("ase", "native"), default="ase",
                        help="Full numerical Hessian: separate ASE gradients or one unbiased native xTB process")
    args = parser.parse_args()
    if not args.xtb or not 0 < args.search_seconds <= 7200 or not 0 <= args.wait_seconds <= 10800:
        parser.error("Specify xTB and bounded search/wait durations")
    if (args.resume_endpoints or args.resume_native_endpoints) and not args.resume_routes:
        parser.error("Endpoint continuation requires explicit --resume-routes")
    if args.resume_endpoints and args.resume_native_endpoints:
        parser.error("Select one endpoint source")
    if args.ts_initial_image is not None and not args.resume_routes:
        parser.error("--ts-initial-image requires explicit --resume-routes")
    if args.resume_last_dimer and (not args.resume_routes or args.resume_endpoints or args.resume_native_endpoints):
        parser.error("--resume-last-dimer requires --resume-routes and cannot be combined with endpoint-only continuation")
    if args.resume_neb_steps is not None and (not args.resume_routes or not 1 <= args.resume_neb_steps <= 1200):
        parser.error("--resume-neb-steps requires --resume-routes and a value from 1 to 1200")
    if args.restart_hessian is not None and (not args.resume_last_dimer or len(args.resume_routes or []) != 1):
        parser.error("--restart-hessian requires one explicit --resume-routes and --resume-last-dimer")
    control = json.loads(args.control.read_text(encoding="utf-8"))
    deadline_epoch = float(control["deadline_epoch"])
    requested = list(dict.fromkeys(args.catalysts[:1] if args.pilot else args.catalysts))
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.resume_routes:
        run_resume_campaign(args, deadline_epoch)
        return
    scheduler_id = f"scheduler_{int(time.time())}_{os.getpid()}"
    scheduler_path = args.output_root / f"{scheduler_id}.json"
    remaining, pending, results = set(requested), {}, []
    wait_deadline = min(deadline_epoch, time.time() + args.wait_seconds)
    status = {"pid": os.getpid(), "started_utc": _utc(), "requested": requested,
              "mode": "pilot" if args.pilot else "campaign", "deadline_epoch": deadline_epoch,
              "workers": args.workers, "threads_per_worker": args.threads, "status": "running"}
    atomic_json(scheduler_path, status)
    print(json.dumps(status), flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        while pending or (remaining and time.time() < deadline_epoch):
            available = _load_best(args.best)
            # The wait allowance limits waiting for absent minima, not dispatch
            # of available catalysts queued behind a longer running search.
            if not pending and time.time() >= wait_deadline and not any(identifier in available for identifier in remaining):
                break
            for identifier in requested:
                if len(pending) >= args.workers or time.time() >= deadline_epoch:
                    break
                if identifier in remaining and identifier in available:
                    remaining.remove(identifier)
                    future = pool.submit(run_one, available[identifier], output_root=str(args.output_root),
                        scratch_root=str(args.scratch_root), executable=args.xtb, threads=args.threads,
                        deadline_epoch=deadline_epoch, search_seconds=min(args.search_seconds, 180) if args.pilot else args.search_seconds,
                        pilot=args.pilot, seed=args.seed + requested.index(identifier) * 83,
                        solvent=args.solvent, solvation_state=args.solvation_state,
                        restrained_product_preconditioning=args.restrained_product_preconditioning,
                        refine_unconverged_band=args.refine_unconverged_band,
                        hessian_method=args.hessian_method)
                    pending[future] = identifier
            if pending:
                completed, _ = wait(pending, timeout=2, return_when=FIRST_COMPLETED)
                for future in completed:
                    identifier = pending.pop(future)
                    try:
                        result = future.result()
                    except Exception as error:
                        result = {"catalyst_id": identifier, "status": "worker_failed", "error": repr(error), "accepted": False}
                    results.append(result)
                    print(json.dumps({"utc": _utc(), "completed": result["catalyst_id"], "status": result["status"]}), flush=True)
            elif remaining:
                time.sleep(2)
            status.update(results=results, waiting_for_active_minima=sorted(remaining), active=list(pending.values()), updated_utc=_utc())
            atomic_json(scheduler_path, status)
    status.update(status="finished", finished_utc=_utc(), missing_active_minima=sorted(remaining), results=results)
    atomic_json(scheduler_path, status)
    print(json.dumps(status), flush=True)


if __name__ == "__main__":
    main()
