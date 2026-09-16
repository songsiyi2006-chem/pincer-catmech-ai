"""Audit selected Phase3 minima without modifying their archived evidence."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from pincer_catmech.kinetics.solvation_sensitivity import (
    checked_source, conditional_association, recalculate_g)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(repo, output):
    if output.exists():
        raise FileExistsError("Use a fresh output; archived results are immutable")
    summary_path = repo / "data/phase3/solvation/summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    sources = {summary_path.relative_to(repo).as_posix(): sha(summary_path)}

    def verify(record):
        path = checked_source(repo, record)
        sources[path.relative_to(repo).as_posix()] = sha(path)
        return path

    def thermo_index(rows):
        indexed = {}
        for row in rows:
            verify({"path": row["source"], "sha256": row["source_sha256"]})
            t = row["temperature"]
            if t in indexed:
                raise ValueError("Duplicate temperature")
            indexed[t] = row
        return indexed

    refs = {}
    for name, record in summary["references"].items():
        for key in ("source_minimum_result", "source_geometry", "source_spectrum"):
            if key not in record:
                raise ValueError(f"Missing reference source: {name}/{key}")
            verify(record[key])
        refs[name] = thermo_index(record["thermochemistry"])
    baseline = {(r["name"], r["temperature_K"]): r["delta_G_qRRHO_1M_kcal_mol"]
                for r in summary["association_thermochemistry"]}
    rows, baseline_errors = [], []
    for cluster in summary["selected"]:
        if (cluster["thermochemistry_status"] != "verified_minimum"
                or cluster["thermochemical_geometry_connectivity_retained"] is not True):
            raise ValueError("Uncertified selected cluster")
        result_path = verify(cluster["stationarity_result"])
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result["accepted"] is not True or result["status"] != "verified_minimum":
            raise ValueError("Minimum certificate rejected")
        geometry = result_path.parent / "accepted_geometry.xyz"
        verify({"path": str(geometry), "sha256": result["accepted_geometry_sha256"]})
        index = thermo_index(result["thermochemistry"])
        if not set(index) <= set(refs["hemiaminal"]) or not set(index) <= set(refs["tert_butanol"]):
            raise ValueError("A cluster temperature has no matching reference; interpolation is forbidden")
        n = cluster["n_solvent"]
        for t, therm in sorted(index.items()):
            for cutoff in (50.0, 100.0, 150.0):
                dg = (recalculate_g(therm, cutoff)
                      - recalculate_g(refs["hemiaminal"][t], cutoff)
                      - n * recalculate_g(refs["tert_butanol"][t], cutoff))
                if cutoff == 100:
                    baseline_errors.append(abs(dg - baseline[(cluster["name"], t)]))
                for activity in (0.001, 0.01, 0.05, 0.1, 1.0):
                    rows.append(dict(name=cluster["name"], n_alcohol=n,
                        temperature_K=t, cutoff_cm1=cutoff, alcohol_activity=activity,
                        delta_G_standard_1M_kcal_mol=dg,
                        **conditional_association(dg, n, t, activity)))
    if len(baseline_errors) != len(baseline) or max(baseline_errors) > 1e-8:
        raise ValueError("Prior baseline was not exactly reconstructed")
    output.mkdir(parents=True)
    with (output / "sensitivity.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    stats = []
    for n in (1, 2, 3):
        subset = [r for r in rows if r["n_alcohol"] == n and r["temperature_K"] == 383.15]
        dgs = [r["delta_G_standard_1M_kcal_mol"] for r in subset]
        stats.append(dict(n_alcohol=n, temperature_K=383.15,
            cutoff_min_kcal_mol=min(dgs), cutoff_max_kcal_mol=max(dgs),
            cutoff_range_kcal_mol=max(dgs)-min(dgs)))
    manifest = dict(schema="pincer_solvation_sensitivity_v1",
        created_utc=datetime.now(timezone.utc).isoformat(), rows=len(rows),
        new_quantum_jobs=0, baseline_rows_checked=len(baseline_errors),
        baseline_max_error_kcal_mol=max(baseline_errors), baseline_temperature_K=383.15,
        cutoff_sensitivity_at_baseline=stats, source_sha256=sources,
        implementation_sha256={p.relative_to(REPO).as_posix(): sha(p) for p in (
            Path(__file__), REPO / "src/pincer_catmech/kinetics/solvation_sensitivity.py",
            REPO / "src/pincer_catmech/kinetics/free_energy.py")},
        csv_sha256=sha(output / "sensitivity.csv"),
        scientific_status="conditional_thermochemical_sensitivity_only",
        assumptions=["One selected certified minimum per cluster size; incomplete conformer ensemble",
            "ALPB electronic solvent contribution held fixed across temperature",
            "Neutral tBuOH activity relative to 1 M is an independent hypothetical input",
            "Alcohol activity is not tBuOK equivalents or a measured solution activity",
            "Association equilibrium only; not catalytic rates, bulk populations or a TS barrier"],
        physical_kinetics_validated=False, bulk_speciation_validated=False)
    (output / "summary.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n",
                                          encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in manifest.items() if k not in {"source_sha256", "implementation_sha256"}}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPO / "data/phase4/solvation_sensitivity_001")
    args = parser.parse_args()
    run(REPO, args.output.resolve())
