"""Collect actual ALPB-toluene Hessians for up to 72 retained catalyst states.

Run with the campaign interpreter; one worker/one CPU thread is intentional.
The generator's best-conformer records are read-only inputs and never changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import time


def audit_existing_results(output_base: Path, records_path: Path):
    """Join minimum certificates to exact historical input identities and hashes."""
    from ase.io import read
    from pincer_catmech.quantum.identity import catalyst_identity_screen
    from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now
    history = {}
    history_bytes = records_path.read_bytes() if records_path.exists() else None
    if history_bytes is not None:
        for line in history_bytes.decode("utf-8").splitlines():
            if not line.strip():
                continue
            item = json.loads(line)
            directory = item.get("native_result", {}).get("output_directory")
            if directory is not None:
                history[os.path.normcase(str(Path(directory).resolve()))] = item["metadata"]
    rows = []
    for path in sorted(output_base.glob("*__*/result.json")):
        result = json.loads(path.read_text(encoding="utf-8"))
        row = {"result_file": str(path), "result_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
               "index_key": result.get("source_metadata", {}).get("index_key"),
               "mathematical_status": result["status"], "mathematical_minimum_accepted": result.get("accepted", False),
               "eligible_for_intended_catalyst_thermochemistry": False}
        if not result.get("finished_utc"):
            row["audit_status"] = "calculation_in_progress"
            rows.append(row)
            continue
        if not result.get("accepted") or result["status"] != "verified_minimum":
            row["audit_status"] = "mathematical_minimum_not_verified"
            rows.append(row)
            continue
        try:
            source = result.get("source_metadata", {}).get("native_geometry")
            metadata = result.get("chemical_identity_metadata") or history[os.path.normcase(str(Path(source).resolve().parent))]
            row["state"] = metadata["state"]
            row["identity_metadata_sha256"] = hashlib.sha256(json.dumps(metadata, sort_keys=True).encode("utf-8")).hexdigest()
            geometry = path.parent / "accepted_geometry.xyz"
            digest = hashlib.sha256(geometry.read_bytes()).hexdigest()
            row.update(accepted_geometry_sha256=digest,
                       geometry_hash_matches_record=digest == result["accepted_geometry_sha256"])
            row["chemical_identity"] = catalyst_identity_screen(read(geometry), metadata)
            source_checks = []
            for record in result["thermochemistry"]:
                spectrum = Path(record["source"])
                source_checks.append(spectrum.exists() and hashlib.sha256(spectrum.read_bytes()).hexdigest() == record["source_sha256"])
            row["thermochemistry_source_hashes_valid"] = bool(source_checks) and all(source_checks)
            row["eligible_for_intended_catalyst_thermochemistry"] = bool(
                result.get("accepted") and result["status"] == "verified_minimum" and
                row["geometry_hash_matches_record"] and row["thermochemistry_source_hashes_valid"] and row["chemical_identity"]["retained"])
            row["audit_status"] = "eligible" if row["eligible_for_intended_catalyst_thermochemistry"] else "identity_or_provenance_rejected"
        except Exception as error:
            row["audit_error"] = f"{type(error).__name__}: {error}"
            row["audit_status"] = "audit_failed"
        rows.append(row)
    audit = {"created_utc": utc_now(), "scope": "Intended catalyst identity is a separate necessary screen from mathematical stationarity; no calculated values are altered.",
             "history_source": str(records_path),
             "history_snapshot_sha256": hashlib.sha256(history_bytes).hexdigest() if history_bytes is not None else None,
             "history_snapshot_bytes": len(history_bytes) if history_bytes is not None else None,
             "rows": rows, "eligible_count": sum(row["eligible_for_intended_catalyst_thermochemistry"] for row in rows)}
    atomic_json(output_base / "identity_audit.json", audit)
    print(json.dumps({"stage": "identity_audit", "audited": len(rows), "eligible": audit["eligible_count"],
                      "in_progress": [row["index_key"] for row in rows if row["audit_status"] == "calculation_in_progress"],
                      "completed_ineligible": [row["index_key"] for row in rows if not row["eligible_for_intended_catalyst_thermochemistry"]
                                               and row["audit_status"] != "calculation_in_progress"]}), flush=True)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best-index", type=Path, default=Path("data/datasets/best_conformers.json"))
    parser.add_argument("--control", type=Path, default=Path("data/campaign/control.json"))
    parser.add_argument("--output-base", type=Path, default=Path("data/campaign/thermochemistry"))
    parser.add_argument("--scratch-base", type=Path, default=Path("../hessian-runs"))
    parser.add_argument("--max-items", type=int, default=24)
    parser.add_argument("--states", nargs="+", choices=("active", "hydrogenated", "protonated_reference"), default=["active"])
    parser.add_argument("--max-repairs", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--catalyst", default=None)
    parser.add_argument("--catalysts", nargs="+", default=None)
    parser.add_argument("--exclude-catalysts", nargs="*", default=[])
    parser.add_argument("--audit-existing", action="store_true")
    parser.add_argument("--records", type=Path, default=Path("data/campaign_alpb_toluene/records.jsonl"))
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.max_items <= 72:
        parser.error("--max-items must be from 1 through 72")
    os.environ.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    from ase.io import read
    from pincer_catmech.quantum.thermochemistry import characterize_stationary_point
    from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now
    if args.audit_existing:
        audit_existing_results(args.output_base, args.records)
        return

    control = json.loads(args.control.read_text(encoding="utf-8"))
    deadline = float(control["deadline_epoch"])
    args.output_base.mkdir(parents=True, exist_ok=True)
    summary = {"status": "running", "started_utc": utc_now(), "worker_count": 1, "threads": 1,
               "solvent": "toluene", "solvation_state": "gsolv", "states": args.states,
               "excluded_catalysts": args.exclude_catalysts,
               "items": [], "deadline_epoch": deadline}
    seen = set()
    while len(seen) < args.max_items and time.time() < deadline:
        index = json.loads(args.best_index.read_text(encoding="utf-8")) if args.best_index.exists() else {}
        candidates = []
        for key, item in index.items():
            identifier = item.get("catalyst_id", "")
            if (item.get("state") not in args.states or not item.get("coordination_retained")
                    or not item.get("native_result", {}).get("converged")
                    or item.get("solvent") != "toluene" or item.get("solvation_state") != "gsolv"
                    or key in seen or identifier in args.exclude_catalysts
                    or (args.catalysts is not None and identifier not in args.catalysts)
                    or (args.catalyst is not None and identifier != args.catalyst)):
                continue
            if not re.fullmatch(r"[A-Za-z0-9_]+", identifier):
                raise ValueError("unexpected catalyst identifier; refusing unsafe output path")
            candidates.append((key, item))
        if not candidates:
            if not args.watch:
                break
            atomic_json(args.output_base / "collector.json", summary)
            time.sleep(min(15, max(0, deadline - time.time())))
            continue
        # Smaller molecules first gives prompt format/force/Hessian feedback.
        candidates.sort(key=lambda pair: (pair[1]["native_result"].get("atoms_count", 10**9), pair[0]))
        for key, item in candidates:
            if len(seen) >= args.max_items or time.time() >= deadline:
                break
            seen.add(key)
            identifier = item["catalyst_id"] + "__" + item["state"]
            output = args.output_base / identifier
            if (output / "result.json").exists():
                result = json.loads((output / "result.json").read_text(encoding="utf-8"))
            else:
                native = item["native_result"]
                geometry = Path(native["output_directory"]) / "xtbopt.xyz"
                print(json.dumps({"utc": utc_now(), "stage": "native_hessian_start", "catalyst": identifier,
                                  "source_geometry": str(geometry), "threads": 1}), flush=True)
                source_hash = hashlib.sha256(geometry.read_bytes()).hexdigest()
                expected_hash = item.get("optimized_geometry_sha256")
                if expected_hash and expected_hash != source_hash:
                    raise ValueError(f"Source geometry changed after indexing: {identifier}")
                atoms = read(geometry)
                metadata = item["metadata"]
                result = characterize_stationary_point(atoms, output, args.scratch_base / identifier,
                    charge=int(native["charge"]), unpaired=int(native["unpaired_electrons"]),
                    threads=1, max_repairs=args.max_repairs, deadline_epoch=deadline,
                    identity_metadata=metadata,
                    solvent="toluene", solvation_state="gsolv",
                    coordination_indices={"metal_index": metadata["metal_index"], "donor_indices": metadata["donor_indices"]},
                    source_metadata={"catalyst_id": item["catalyst_id"], "state": item["state"], "index_key": key,
                                     "conformer": item.get("conformer"), "geometry_sha256": source_hash,
                                     "source_index": str(args.best_index.resolve()), "native_geometry": str(geometry),
                                     "spin_provenance": metadata.get("spin_provenance"),
                                     "formula": metadata.get("formula")})
            record = {"key": key, "status": result["status"], "accepted": result["accepted"],
                      "result_file": str(output / "result.json"), "wall_seconds": result.get("wall_seconds")}
            summary["items"].append(record)
            atomic_json(args.output_base / "collector.json", summary)
            print(json.dumps({"utc": utc_now(), "stage": "native_hessian_complete", **record}), flush=True)
    summary.update(status="deadline_reached" if time.time() >= deadline else "completed_available_inputs",
                   finished_utc=utc_now(), accepted_minima=sum(item["accepted"] for item in summary["items"]))
    atomic_json(args.output_base / "collector.json", summary)
    print(json.dumps(summary), flush=True)
    audit_existing_results(args.output_base, args.records)


if __name__ == "__main__":
    main()
