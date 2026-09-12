"""Rebuild a compact audit of actual condensation searches and endpoint minima.

Every endpoint certification is checked against its saved full Hessian before
writing a standalone record. Failed reaction paths remain failed.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/"src"))
import numpy as np
from ase.io import read
from pincer_catmech.quantum.xtb_backend import atomic_json


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    directory = REPO/"data"/"condensation_search"
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed_bounded_searches", "method": "GFN2-xTB",
        "solvent": "ALPB toluene", "solvation_state": "gsolv",
        "energy_origin": "Actual native xTB energies, gradients and Hessians only",
        "scope": "Explicit cluster models; no potassium, no total-base kinetic law, no experimental validation",
        "attempts": [], "verified_endpoint_minimum_records": [],
        "synthetic_research_observations": 0, "accepted_transition_states": 0,
        "complete_kinetic_network": False, "barriers_and_TOFs_available": False}
    for result_path in sorted(directory.glob("*/result.json")):
        record = json.loads(result_path.read_text(encoding="utf-8"))
        folder = result_path.parent
        row = {"id": folder.name, "status": record["status"], "accepted": record["accepted"],
            "started_utc": record["started_utc"], "finished_utc": record.get("finished_utc"),
            "elapsed_seconds": record.get("elapsed_seconds"), "charge": record["charge"],
            "model": record["model"], "result_path": str(result_path), "result_sha256": sha(result_path),
            "failure_reasons": record["failure_reasons"], "endpoint_identities": [
                {"index": endpoint["index"], **endpoint["identity"]} for endpoint in record["endpoints"]]}
        logfile = folder/"neb.log"
        if logfile.exists():
            values = [line.split() for line in logfile.read_text().splitlines() if line.startswith("FIRE:")]
            if values:
                eligible = [value for value in values if int(value[1]) >= 31]
                best = min(eligible, key=lambda value: float(value[-1])) if eligible else None
                row.update(neb_steps=int(values[-1][1]), best_post_warmup_step=int(best[1]) if best else None,
                    minimum_post_warmup_projected_fmax_eV_A=float(best[-1]) if best else None,
                    final_projected_fmax_eV_A=float(values[-1][-1]), n_total_images=9)
            if (folder/"neb.traj").exists():
                row.update(trajectory_path=str(folder/"neb.traj"), trajectory_sha256=sha(folder/"neb.traj"))
        summary["attempts"].append(row)
        summary["accepted_transition_states"] += int(record["accepted"])
        for endpoint in record["endpoints"]:
            if endpoint.get("stationarity_status") != "verified_minimum":
                continue
            index = endpoint["index"]
            hessian_path, geometry_path = folder/f"endpoint_{index}_hessian.npz", folder/f"endpoint_{index}.xyz"
            frequencies = np.asarray(endpoint["frequencies_cm1"], dtype=float)
            atoms = read(geometry_path)
            with np.load(hessian_path) as saved:
                if not np.array_equal(atoms.numbers, saved["atomic_numbers"]) or not np.allclose(atoms.positions, saved["positions"], rtol=0, atol=1e-6):
                    raise ValueError("Endpoint geometry differs from its Hessian artifact")
                if not np.allclose(frequencies, saved["frequencies_cm1"], rtol=0, atol=1e-8):
                    raise ValueError("Endpoint reported spectrum differs from its Hessian")
                if len(frequencies) != 3*len(atoms)-int(saved["external_rank"]) or np.any(frequencies <= 0):
                    raise ValueError("Endpoint lacks a complete all-positive internal spectrum")
            if sha(hessian_path) != endpoint["hessian_sha256"] or endpoint["hessian_antisymmetry_relative"] > 0.05:
                raise ValueError("Endpoint Hessian hash or numerical-consistency criterion failed")
            standalone = {"status": "verified_minimum", "accepted": True,
                "scope": "Endpoint minimum only; this does not certify the containing reaction path",
                "method": record["method"], "solvent": record["solvent"], "solvation_state": record["solvation_state"],
                "model": record["model"], "charge": record["charge"], "unpaired_electrons": 0,
                "source_path": str(hessian_path), "source_sha256": sha(hessian_path),
                "source_geometry": str(geometry_path), "source_geometry_sha256": sha(geometry_path),
                "parent_result": str(result_path), "parent_result_sha256": sha(result_path),
                "parent_path_status": record["status"], "frequencies_cm1": frequencies.tolist(),
                "minimum_frequency_cm1": float(frequencies.min()), "internal_modes": len(frequencies),
                "thermochemistry": endpoint["thermochemistry"]}
            target = folder/f"endpoint_{index}_verified_minimum.json"
            atomic_json(target, standalone)
            summary["verified_endpoint_minimum_records"].append({"path": str(target),
                "sha256": sha(target), "model": record["model"], "minimum_frequency_cm1": float(frequencies.min()),
                "internal_modes": len(frequencies), "formula": atoms.get_chemical_formula(),
                "charge": record["charge"]})
    summary["completed_search_attempts"] = len(summary["attempts"])
    summary["unique_verified_endpoint_chemical_models"] = len({row["model"] for row in summary["verified_endpoint_minimum_records"]})
    summary["full_endpoint_hessian_records"] = len(summary["verified_endpoint_minimum_records"])
    if any(row["status"] == "running" for row in summary["attempts"]):
        summary["status"] = "running"
    atomic_json(directory/"summary.json", summary)
    print(json.dumps({key: summary[key] for key in ("status", "completed_search_attempts", "accepted_transition_states",
        "full_endpoint_hessian_records", "unique_verified_endpoint_chemical_models")}))


if __name__ == "__main__":
    main()
