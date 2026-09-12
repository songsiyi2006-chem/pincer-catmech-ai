"""One full native Hessian diagnostic of a completed search's actual last dimer frame.

This operator utility never certifies a TS or calculates any G/barrier.
Usage: diagnose_ts_candidate.py COMPLETED_ROUTE [SECONDS [INTERNAL_IMAGE]].
The optional image selects an actual frame from the final complete nine-image
band; otherwise the last actual dimer frame is used. All calculations retain
the campaign control file's deadline and require complete chemical identity.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import time

import numpy as np
from ase.io import write
from ase.data import covalent_radii
from ase.io.trajectory import Trajectory

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from run_neb_campaign import make_native_hessian_provider, ligand_identity_screen
from pincer_catmech.kinetics.neb_ts_search import transfer_mode_overlap
from pincer_catmech.quantum.xtb_backend import XTBCalculator, atomic_json
from pincer_catmech.generators import CatalystSpec, generate_catalyst

source = Path(sys.argv[1]).resolve()
task_source = json.loads((source / "task.json").read_text(encoding="utf-8"))
path_result = json.loads((source / "path/result.json").read_text(encoding="utf-8"))
assembly = json.loads((source / "assembly.json").read_text(encoding="utf-8"))
if task_source.get("solvent") != "toluene" or task_source.get("solvation_state") != "gsolv":
    raise ValueError("Diagnostic requires an actual ALPB toluene/gsolv source")
candidate = None
selected_image = int(sys.argv[3]) if len(sys.argv) > 3 else None
if selected_image is not None and not 1 <= selected_image <= 7:
    raise ValueError("Band diagnostic must select an explicit internal image 1 through 7")
for trajectory_source in ([source / "path/neb.traj"] if selected_image is not None else
                          sorted((source / "path").glob("attempt_*/dimer.traj"), reverse=True)):
    with Trajectory(trajectory_source, "r") as trajectory:
        if len(trajectory):
            if selected_image is not None and len(trajectory) // 9 < 1:
                raise ValueError("Source has no actual complete nine-image band")
            frame = (len(trajectory) // 9 - 1) * 9 + selected_image if selected_image is not None else len(trajectory) - 1
            candidate = trajectory[frame]
            break
if candidate is None:
    raise ValueError("No actual saved dimer candidate is available")
candidate.info = dict(assembly["electronic_state"])
metadata = assembly.get("catalyst_metadata")
if metadata is None:
    metal, remainder = task_source["catalyst_id"].split("_", 1)
    backbone, substituent = remainder.rsplit("_", 1)
    metadata = generate_catalyst(CatalystSpec(metal, backbone, substituent), state="active").metadata()
from rdkit import Chem
from pincer_catmech.quantum.identity import expected_ligand_graph
bonds, _ = expected_ligand_graph(*(metadata[key] for key in ("metal", "backbone", "substituent", "state")))
required = {tuple(sorted(pair)) for pair in bonds} | {tuple(sorted(pair)) for pair in metadata["carbonyl_indices"]}
substrate = Chem.AddHs(Chem.MolFromSmiles("[OH:1][CH2:2]c1ccccc1"))
offset = assembly["catalyst_atom_count"]
required |= {tuple(sorted((offset + bond.GetBeginAtomIdx(), offset + bond.GetEndAtomIdx()))) for bond in substrate.GetBonds()}
transfer = assembly["transfer_indices"]
variable = {tuple(sorted((transfer[channel], transfer[channel + "_donor"]))) for channel in ("proton", "hydride")}
required -= variable
variable.add(tuple(sorted((transfer["proton"], transfer["proton_acceptor"]))))
observed = {(j, i) for i in range(1, len(candidate)) for j in range(1, i)
            if candidate.get_distance(i, j) < 1.28 * (covalent_radii[candidate[i].number] + covalent_radii[candidate[j].number])}
full_graph = {"retained": not (required - observed or observed - required - variable),
              "lost_nonreacting_bonds": sorted(required - observed),
              "unexpected_nonmetal_bonds": sorted(observed - required - variable),
              "allowed_variable_transfer_pairs": sorted(variable)}
control = json.loads((REPO / "data/campaign/control.json").read_text(encoding="utf-8"))
duration = float(sys.argv[2]) if len(sys.argv) > 2 else 900
if not np.isfinite(duration) or not 0 < duration <= 1800:
    raise ValueError("Diagnostic duration must be positive and at most 1800 seconds")
deadline = min(control["deadline_epoch"], time.time() + duration)
label = f"{task_source['catalyst_id']}_diagnostic_{int(time.time())}_{os.getpid()}"
directory = REPO / "data/campaign/neb" / label
directory.mkdir(exist_ok=False)
write(directory / "candidate.xyz", candidate, format="xyz")
record = {"status": "started", "accepted": False, "scope": "Full nonstationary-candidate curvature diagnostic; no TS certificate, thermochemistry, or barrier is emitted",
    "catalyst_id": task_source["catalyst_id"], "source_route": str(source), "source_result_status": path_result["status"],
    "source_trajectory": str(trajectory_source), "source_trajectory_sha256": hashlib.sha256(trajectory_source.read_bytes()).hexdigest(),
    "source_frame": frame, "candidate_xyz_sha256": hashlib.sha256((directory / "candidate.xyz").read_bytes()).hexdigest(),
    "selected_band_image": selected_image,
    "selection_reason": "Explicit requested actual internal band image; diagnostic significance must be assessed from the archived mapped distances" if selected_image is not None else "Actual last saved dimer candidate",
    "full_nonreacting_graph": full_graph,
    "transfer_distances_A": {channel: {"donor_H": candidate.get_distance(transfer[channel], transfer[channel + "_donor"]),
        "acceptor_H": candidate.get_distance(transfer[channel], transfer[channel + "_acceptor"])} for channel in ("proton", "hydride")},
    "started_utc": datetime.now(timezone.utc).isoformat(), "pid": os.getpid(), "threads": 2,
    "chemical_identity": ligand_identity_screen(candidate, metadata), "solvent": "toluene", "solvation_state": "gsolv"}
atomic_json(directory / "diagnostic.json", record)
print(json.dumps({"directory": str(directory), "status": "started", "pid": os.getpid()}), flush=True)
try:
    if not record["chemical_identity"]["retained"] or not full_graph["retained"]:
        raise ValueError("Fixed candidate fails full catalyst/nonreacting-substrate identity; diagnostic calculation skipped")
    candidate.calc = XTBCalculator(REPO.parent / "neb-scratch" / label / "gradient",
        charge=candidate.info["charge"], unpaired=candidate.info["unpaired"], threads=2,
        executable=os.environ.get("PINCER_XTB"), deadline=deadline, solvent="toluene", solvation_state="gsolv",
        trace_path=directory / "fresh_gradient.jsonl.gz")
    record["true_fmax_eV_A"] = float(np.linalg.norm(candidate.get_forces(), axis=1).max())
    provider = make_native_hessian_provider(REPO.parent / "neb-scratch", label,
        charge=candidate.info["charge"], unpaired=candidate.info["unpaired"], threads=2,
        executable=os.environ.get("PINCER_XTB"), deadline=deadline, solvent="toluene", solvation_state="gsolv")
    spectrum = provider(candidate, directory)
    spectrum.save(directory / "diagnostic_full_hessian.npz", candidate)
    negative = np.flatnonzero(spectrum.frequencies_cm1 < 0)
    record.update(status="diagnostic_complete", curvature_frequencies_cm1=spectrum.frequencies_cm1.tolist(),
        negative_mode_count=len(negative), external_rank=spectrum.external_rank,
        hessian_antisymmetry_relative=spectrum.antisymmetry_relative,
        negative_mode_transfer_coordinates=[{"frequency_cm1": float(spectrum.frequencies_cm1[i]),
            **transfer_mode_overlap(candidate, spectrum.cartesian_modes[i], assembly["transfer_indices"])} for i in negative],
        hessian_npz_sha256=hashlib.sha256((directory / "diagnostic_full_hessian.npz").read_bytes()).hexdigest())
except Exception as error:
    record.update(status="diagnostic_failed", error=f"{type(error).__name__}: {error}")
record["finished_utc"] = datetime.now(timezone.utc).isoformat()
atomic_json(directory / "diagnostic.json", record)
print(json.dumps({key: value for key, value in record.items() if key not in ("curvature_frequencies_cm1", "chemical_identity")}), flush=True)
