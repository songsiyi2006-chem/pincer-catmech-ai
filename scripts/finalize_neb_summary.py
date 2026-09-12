"""Summarize actual completed NEB evidence without manufacturing energy labels."""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "data/campaign/neb"
sys.path.insert(0, str(REPO / "src"))
from pincer_catmech.quantum.xtb_backend import atomic_json


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else {}


def successful_gradients(folder):
    rows = []
    for path in sorted(set(folder.glob("evaluations_*.jsonl.gz")) | set(folder.glob("fresh_gradient.jsonl.gz"))):
        count, errors = 0, []
        try:
            with gzip.open(path, "rt", encoding="utf-8") as stream:
                for number, line in enumerate(stream, 1):
                    try:
                        record = json.loads(line)
                        if record.get("forces_eV_A") is not None and record.get("energy_eV") is not None:
                            count += 1
                    except Exception as error:
                        errors.append(f"line {number}: {type(error).__name__}: {error}")
        except Exception as error:
            errors.append(f"archive: {type(error).__name__}: {error}")
        rows.append({"file": path.name, "sha256": digest(path), "successful_gradient_records": count, "read_errors": errors})
    return rows


def stage_history(folder):
    source = folder / "path/progress.jsonl"
    stages, errors = [], []
    for number, line in enumerate(source.read_text(encoding="utf-8").splitlines() if source.exists() else [], 1):
        try:
            stage = json.loads(line).get("stage")
            if stage not in stages:
                stages.append(stage)
        except Exception as error:
            errors.append(f"line {number}: {type(error).__name__}: {error}")
    return {"stage_history": stages, "progress_parse_errors": errors}


def accepted_identity(folder):
    """Check fixed accepted coordinates while allowing only mapped transfer bonds to vary."""
    from ase.io import read
    from ase.data import covalent_radii
    from rdkit import Chem
    from pincer_catmech.quantum.identity import catalyst_identity_screen, expected_ligand_graph
    source = folder / "path/accepted_ts.xyz"
    if not source.exists():
        return {"eligible": False, "reason": "Engine acceptance has no fixed accepted_ts.xyz"}
    assembly = load(folder / "assembly.json")
    metadata = assembly.get("catalyst_metadata") or load(REPO / "data/datasets/best_conformers.json")[assembly["catalyst_id"] + "|active"]["metadata"]
    atoms = read(source)
    identity = catalyst_identity_screen(atoms, metadata)
    bonds, _ = expected_ligand_graph(*(metadata[key] for key in ("metal", "backbone", "substituent", "state")))
    expected = {tuple(sorted(pair)) for pair in bonds} | {tuple(sorted(pair)) for pair in metadata["carbonyl_indices"]}
    substrate = Chem.AddHs(Chem.MolFromSmiles("[OH:1][CH2:2]c1ccccc1"))
    offset = assembly["catalyst_atom_count"]
    expected |= {tuple(sorted((offset + bond.GetBeginAtomIdx(), offset + bond.GetEndAtomIdx()))) for bond in substrate.GetBonds()}
    transfer = assembly["transfer_indices"]
    variable = {tuple(sorted((transfer[channel], transfer[channel + "_donor"]))) for channel in ("proton", "hydride")}
    expected -= variable
    variable.add(tuple(sorted((transfer["proton"], transfer["proton_acceptor"]))))
    observed = {(j, i) for i in range(1, len(atoms)) for j in range(1, i)
                if atoms.get_distance(i, j) < 1.28 * (covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number])}
    lost, new = expected - observed, observed - expected - variable
    return {"eligible": bool(identity["retained"] and not lost and not new),
            "accepted_geometry_sha256": digest(source), "catalyst_identity": identity,
            "required_nonmetal_bonds_lost": sorted(lost), "unexpected_nonmetal_bonds": sorted(new),
            "allowed_variable_transfer_pairs": sorted(variable),
            "scope": "Original catalyst and all nonreacting substrate covalent connectivity; only the mapped O-H/C-H/N-H transfer pairs may change at the saddle"}


query = "Get-CimInstance Win32_Process | Where-Object { $_.Name -match '^python(w)?\\.exe$' -and $_.CommandLine -match '(run_neb_campaign|diagnose_ts_candidate)\\.py' } | Select-Object ProcessId,ParentProcessId,CommandLine | ConvertTo-Json -Compress"
snapshot = subprocess.run(["powershell", "-NoProfile", "-Command", query], check=True, text=True, capture_output=True).stdout.strip()
active = json.loads(snapshot) if snapshot else []
if active:
    raise RuntimeError(f"Refusing final summary while tracked NEB writers still run: {active}")

excluded_file = ROOT / "excluded_route_identities.json"
exclusions = {row["route"]: row for row in load(excluded_file).get("records", [])}
identifiers = [f"{metal}_macho_pnp_{substituent}" for metal in ("Ru", "Mn", "Co") for substituent in ("Ph", "iPr")]
groups = {identifier: {"catalyst_id": identifier, "routes": [], "diagnostics": []} for identifier in identifiers}
for folder in sorted(ROOT.iterdir()):
    if not folder.is_dir():
        continue
    task, diagnostic = load(folder / "task.json"), load(folder / "diagnostic.json")
    record = task or diagnostic
    identifier = record.get("catalyst_id")
    if identifier not in groups:
        continue
    traces = successful_gradients(folder)
    row = {"directory": folder.name, "status": record.get("status"), "accepted": bool(record.get("accepted", False)),
           "started_utc": record.get("started_utc"), "finished_utc": record.get("finished_utc"),
           "successful_archived_gradient_records": sum(item["successful_gradient_records"] for item in traces),
           "gradient_archives": traces, "solvent": record.get("solvent"), "solvation_state": record.get("solvation_state"),
           "source_route": record.get("source_route"), "identity_exclusion": exclusions.get(folder.name),
           "source_geometry_sha256": record.get("source_geometry_sha256"),
           "source_trajectory_sha256": record.get("source_trajectory_sha256"),
           "source_endpoint_sha256": record.get("source_endpoint_sha256"),
           "task_or_diagnostic_sha256": digest(folder / ("task.json" if task else "diagnostic.json"))}
    if diagnostic:
        provider = load(folder / "native_hessian_provider.json")
        row.update(scope=diagnostic.get("scope"), true_fmax_eV_A=diagnostic.get("true_fmax_eV_A"),
                   negative_mode_count=diagnostic.get("negative_mode_count"),
                   negative_mode_transfer_coordinates=diagnostic.get("negative_mode_transfer_coordinates"),
                   chemical_identity=diagnostic.get("chemical_identity"),
                   full_nonreacting_graph=diagnostic.get("full_nonreacting_graph"),
                   source_frame=diagnostic.get("source_frame"), selected_band_image=diagnostic.get("selected_band_image"),
                   selection_reason=diagnostic.get("selection_reason"), transfer_distances_A=diagnostic.get("transfer_distances_A"),
                   full_matrix_shape=provider.get("full_matrix_shape"),
                   native_gradient_energy_difference_eV=provider.get("energy_difference_eV"),
                   native_hessian_provider_sha256=digest(folder / "native_hessian_provider.json") if provider else None,
                   hessian_npz_sha256=diagnostic.get("hessian_npz_sha256"), error=diagnostic.get("error"))
        groups[identifier]["diagnostics"].append(row)
        continue
    result, protocol = load(folder / "path/result.json"), load(folder / "path/protocol.json")
    final_identity = accepted_identity(folder) if result.get("accepted") else None
    row["engine_accepted"] = bool(result.get("accepted", False))
    row["accepted_ts_identity"] = final_identity
    row["source_not_explicitly_excluded"] = not bool(row["identity_exclusion"])
    row["intended_identity_eligible"] = (False if row["identity_exclusion"] else
        final_identity["eligible"] if final_identity is not None else None)
    row["accepted"] = row["engine_accepted"] and row["intended_identity_eligible"] is True
    row.update(path_status=result.get("status"), endpoint_converged=result.get("endpoint_converged"),
               neb_converged=result.get("neb_converged"), provisional_neb_initial_guess=result.get("provisional_neb_initial_guess"),
               source_identity=task.get("source_ligand_identity"), endpoint_identity=task.get("endpoint_ligand_identity"),
               hessian_mode_restart=task.get("hessian_mode_restart"), attempts=result.get("attempts", []),
               failure_reasons=result.get("failure_reasons", []) + ([task["error"]] if task.get("error") else []),
               full_hessian_files=[str(path.relative_to(folder)) for path in folder.glob("path/attempt_*/full_hessian.npz")],
               **stage_history(folder),
               required_gates={key: protocol.get(key) for key in ("ts_fmax_ev_a", "imaginary_noise_cm1", "imaginary_window_cm1", "minimum_transfer_overlap", "minimum_combined_overlap")},
               forward_electronic_barrier_ev=result.get("forward_electronic_barrier_ev"),
               raw_engine_forward_electronic_barrier_ev=result.get("forward_electronic_barrier_ev"),
               free_energy_barrier_available=False,
               result_sha256=digest(folder / "path/result.json") if result else None)
    if row["intended_identity_eligible"] is False:
        row["failure_reasons"].append("Source or fixed accepted geometry fails the full intended catalyst/substrate identity screen")
        row["forward_electronic_barrier_ev"] = None
    groups[identifier]["routes"].append(row)

finished = datetime.now(timezone.utc).isoformat()
for group in groups.values():
    group["route_count"] = len(group["routes"])
    group["diagnostic_count"] = len(group["diagnostics"])
    group["successful_archived_gradient_records"] = sum(row["successful_archived_gradient_records"] for row in group["routes"] + group["diagnostics"])
    group["accepted_route_count"] = sum(row["accepted"] for row in group["routes"])
    group["identity_excluded_route_count"] = sum(bool(row["identity_exclusion"]) for row in group["routes"])
all_rows = [row for group in groups.values() for key in ("routes", "diagnostics") for row in group[key]]
summary = {"schema_version": 1, "finished_utc": finished,
    "scope": "Actual GFN2-xTB ALPB(toluene)/gsolv six-representative alcohol dehydrogenation searches; software tests excluded",
    "all_tracked_actual_neb_calculation_sessions_ended": True,
    "process_guard": {"checked_utc": finished, "active_python_neb_or_diagnostic_writers": active,
                      "scope": "Tracked run_neb_campaign.py and diagnose_ts_candidate.py sessions; parent finalizer additionally checks all scientific writers"},
    "gradient_count_definition": "Successful JSONL records with both an energy and full forces; excludes native optimization/Hessian internal evaluations and failed SCC attempts; not total xTB calls",
    "successful_archived_gradient_records": sum(row["successful_archived_gradient_records"] for row in all_rows),
    "accepted_path_count": sum(row["accepted"] for group in groups.values() for row in group["routes"]),
    "free_energy_barrier_count": 0,
    "diagnostic_boundary": "Nonstationary diagnostic frequencies and modes do not establish TS identity, a free energy, a barrier, or kinetic labels",
    "identity_exclusion_file": {"file": excluded_file.name, "sha256": digest(excluded_file)},
    "chemical_identity_audits": [{"file": path.name, "sha256": digest(path)} for path in
        (ROOT / "current_reaction_complete_nonmetal_graph_audit.json", ROOT / "co_retry_endpoint_audit.json",
         ROOT / "ru_conformer_identity_index_clarification.json") if path.exists()],
    "groups": list(groups.values())}
atomic_json(ROOT / "FINAL_NEB_SUMMARY.json", summary)
print(json.dumps({key: value for key, value in summary.items() if key != "groups"}))
