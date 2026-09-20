"""Regenerate results strictly from immutable native single-point evidence."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

HARTREE_TO_EV = 27.211386245988


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_pair_protocol(rn, ra):
    if rn["psi4_version"] != ra["psi4_version"]:
        raise ValueError("Different Psi4 versions in charge pair")
    for key in ["method", "basis", "solvent", "geometry_optimization", "frequency_calculation"]:
        if rn["settings"][key] != ra["settings"][key]:
            raise ValueError(f"Incompatible energy-pair setting: {key}")
    if rn["input_sha256"] != ra["input_sha256"]:
        raise ValueError("A fixed-nuclei pair must share an identical geometry hash")
    # Reference must change with spin, and an explicitly recorded SOSCF rescue
    # changes the iterative solver only, not the Hamiltonian/grid or tolerances.
    solver_only = {"reference", "soscf", "soscf_start_convergence", "soscf_max_iter"}
    neutral_options = {k: v for k, v in rn["psi4_options"].items() if k not in solver_only}
    anion_options = {k: v for k, v in ra["psi4_options"].items() if k not in solver_only}
    if neutral_options != anion_options:
        raise ValueError("Incompatible SCF/DF/grid protocol in charge pair")
    if rn["basis_source_sha256"] != ra["basis_source_sha256"]:
        raise ValueError("Different basis source content in charge pair")


def summarize(root: Path):
    rows = []
    validated = {}
    exclusion_path = root / "excluded_scratch_manifest.json"
    previous_exclusions = json.loads(exclusion_path.read_text())["files"] if exclusion_path.exists() else []
    excluded = {r["path"]: r for r in previous_exclusions}
    for folder in sorted((root / "runs").iterdir()):
        if not folder.is_dir():
            continue
        record = json.loads((folder / "run_record.json").read_text())
        for name, digest in record.get("files", {}).items():
            relative = Path(name)
            if "scratch" in relative.parts or relative.name.startswith("psi."):
                continue  # explicitly rebuildable; absent in Git-delivered archive
            if sha256(folder / name) != digest:
                raise ValueError(f"Modified archived output: {folder.name}/{name}")
        row = {
            "run_id": folder.name, "status": record["status"],
            "elapsed_seconds": record.get("elapsed_seconds"),
            "peak_rss_mib": record.get("sampled_peak_process_tree_rss_bytes", 0) / 1024**2,
            "scf_stability_analysis_performed": False,
            "coordinate_source": "Inherited, independently generated GFN2-xTB neutral minimum; not SI coordinates",
            "method": record["config"]["method"], "basis": record["config"]["basis"],
            "charge": record["config"]["charge"],
        }
        for p in folder.rglob("*"):
            if p.is_file() and ("scratch" in p.relative_to(folder).parts or p.name.startswith("psi.")):
                relative_name = str(p.relative_to(root)).replace("\\", "/")
                excluded[relative_name] = {
                    "path": relative_name,
                    "bytes": p.stat().st_size, "sha256": sha256(p),
                    "reason": "Rebuildable native scratch; retained locally, excluded from Git delivery",
                }
        if row["status"] != "completed":
            rows.append(row)
            continue
        result = json.loads((folder / "result.json").read_text())
        if result["status"] != "completed":
            raise ValueError("Runner succeeded but quantum result did not")
        native = (folder / "psi4.out").read_text(encoding="utf-8", errors="replace")
        matches = re.findall(r"@DF-[RU]KS Final Energy:\s*([-+0-9.]+)", native)
        if not matches or "Energy and wave function converged." not in native:
            raise ValueError("Missing native SCF-convergence/final-energy evidence")
        if abs(float(matches[-1]) - result["energy_hartree"]) > 1e-10:
            raise ValueError("Native/result energy mismatch")
        settings = result["settings"]
        settings_file = folder / "settings.json"
        if settings != json.loads(settings_file.read_text()) or settings != record["config"]:
            raise ValueError("Inconsistent result/settings-file/runner protocol fields")
        if sha256(settings_file) != result["settings_sha256"]:
            raise ValueError("Embedded settings hash mismatch")
        if sha256(folder / "executed_psi4_single_point.py") != result["child_script_sha256"]:
            raise ValueError("Embedded child-script hash mismatch")
        expected_reference = "rks" if settings["multiplicity"] == 1 else "uks"
        if result["psi4_options"]["reference"].lower() != expected_reference:
            raise ValueError("Spin reference disagrees with multiplicity")
        for key in ["energy_hartree", "alpha_homo_hartree", "beta_homo_hartree", "s2_ideal"]:
            if not math.isfinite(result[key]):
                raise ValueError(f"Nonfinite result field: {key}")
        if not math.isfinite(float(matches[-1])):
            raise ValueError("Nonfinite native final energy")
        if sha256(folder / "input.xyz") != result["input_sha256"]:
            raise ValueError("Geometry identity mismatch")
        atoms = (folder / "input.xyz").read_text().splitlines()[2:]
        z = {"H": 1, "C": 6, "N": 7, "O": 8}
        electron_count = sum(z[a.split()[0]] for a in atoms if a.strip()) - settings["charge"]
        if electron_count != result["nalpha"] + result["nbeta"]:
            raise ValueError("Electron count mismatch")
        if result["nalpha"] - result["nbeta"] != settings["multiplicity"] - 1:
            raise ValueError("Spin occupation mismatch")
        ideal_s2 = (settings["multiplicity"] - 1) * (settings["multiplicity"] + 1) / 4
        if abs(result["s2_ideal"] - ideal_s2) > 1e-12:
            raise ValueError("Ideal spin expectation disagrees with multiplicity")
        if not math.isfinite(result["s2_orbital_overlap"]) or result["s2_orbital_overlap"] < -1e-8:
            raise ValueError("Invalid determinant spin expectation")
        row.update({key: result[key] for key in [
            "energy_hartree", "nalpha", "nbeta", "nbf", "s2_orbital_overlap",
            "s2_ideal", "alpha_homo_hartree", "beta_homo_hartree", "input_sha256",
        ]})
        row.update(method=settings["method"], basis=settings["basis"], charge=settings["charge"], native_scf_converged=True)
        row["spin_contamination_s2_minus_ideal"] = result["s2_orbital_overlap"] - ideal_s2
        rows.append(row)
        state_key = (settings["method"], settings["charge"])
        if state_key in validated:
            raise ValueError("Ambiguous duplicate successful method/charge state; define an accepted-run registry")
        validated[state_key] = (row, result)
    pairs = []
    for method in ["pbe0", "b3lyp"]:
        if (method, 0) not in validated or (method, -1) not in validated:
            continue
        neutral, rn = validated[(method, 0)]
        anion, ra = validated[(method, -1)]
        validate_pair_protocol(rn, ra)
        dn = (anion["energy_hartree"] - neutral["energy_hartree"]) * HARTREE_TO_EV
        if not math.isfinite(dn):
            raise ValueError("Nonfinite energy difference")
        pairs.append({
            "method": method, "basis": neutral["basis"],
            "neutral_run_id": neutral["run_id"], "anion_run_id": anion["run_id"],
            "delta_E_anion_minus_neutral_eV": dn,
            "electronic_attachment_energy_neutral_minus_anion_eV": -dn,
            "sign_convention": "Positive E(neutral)-E(anion) means the anion is lower in this finite-basis electronic model at the same nuclei; it is not a calibrated reduction potential or proof of a physical bound state",
            "evidence_class": "CALCULATION",
            "anion_alpha_homo_hartree": anion["alpha_homo_hartree"],
            "bound_state_warning": "Finite diffuse basis does not establish a bound anion; positive occupied orbital energy and/or negative attachment energy require extra caution",
            "anion_solver_rescue": ra["settings"].get("soscf_rescue", False),
        })
    external_preflights = 0
    external_path = root / "memory_preflight_rejection.json"
    if external_path.exists():
        external = json.loads(external_path.read_text())
        if external.get("status") != "resource_preflight_not_launched" or external.get("native_psi4_process_launched") is not False:
            raise ValueError("Invalid external preflight event")
        same_event_in_runs = any(r["run_id"] == external.get("requested_job") and r["status"] == "resource_preflight_not_launched" for r in rows)
        external_preflights = 0 if same_event_in_runs else 1
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_class": "CALCULATION", "scope": "molecular gas-phase fixed-nuclei DFT diagnostic",
        "hartree_to_eV": HARTREE_TO_EV, "runs": rows, "energy_pairs": pairs,
        "completed_quantum_runs": sum(r["status"] == "completed" for r in rows),
        "failed_or_terminated_runs": sum(r["status"] in ["failed", "timeout", "memory_limit", "launch_failed", "monitor_error"] for r in rows),
        "preflight_rejections_not_counted_as_quantum_runs": external_preflights + sum(r["status"] == "resource_preflight_not_launched" for r in rows),
        "method_spread_eV": abs(pairs[0]["delta_E_anion_minus_neutral_eV"] - pairs[1]["delta_E_anion_minus_neutral_eV"]) if len(pairs) == 2 else None,
        "limitations": [
            "No DFT geometry optimization, Hessian, stationary-point certification or SCF stability analysis",
            "Single SAD guess per state; orbital-energy signs do not establish bound-state character",
            "One diffuse basis, one grid and one geometry; basis/grid/geometry convergence untested",
            "Two hybrid functionals can share systematic error; agreement is not calibration",
            "No Cu site, periodic interface, solvent, electrode reference, electron chemical potential or electrochemical potential",
            "No xTB electronic energy is mixed into any DFT energy difference",
            "Sampled RSS monitor may miss sub-sampling memory spikes; total sampled peak is reported",
        ],
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    exclusion_path.write_text(json.dumps({"files": sorted(excluded.values(), key=lambda r: r["path"]), "policy": "Retain locally; omit from Git; no scientific label depends on scratch"}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["completed_quantum_runs", "failed_or_terminated_runs", "energy_pairs", "method_spread_eV"]}, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    summarize(args.root)
