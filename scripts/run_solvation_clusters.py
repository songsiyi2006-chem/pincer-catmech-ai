"""Sequential, evidence-preserving GFN2-xTB/ALPB(toluene) microsolvation campaign."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
import zipfile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from ase.io import read, write
from ase.units import Hartree
import numpy as np

from pincer_catmech.generators.solvation_clusters import (
    HEMIAMINAL_SMILES, TERT_BUTANOL_SMILES, EV_TO_KCAL_MOL,
    association_energy, build_solvation_cluster, connectivity_check,
    graph_bonds, hydrogen_bond_metrics, many_body_nonadditivity, molecular_graph,
)
from pincer_catmech.quantum.thermochemistry import characterize_stationary_point, thermochemistry_from_modes
from pincer_catmech.quantum.xtb_backend import (
    _execute, atomic_json, executable_path, optimize_geometry, parse_energy, utc_now,
)

PROTOCOL = {"method": "GFN2-xTB", "solvent": "toluene", "solvation_state": "gsolv",
            "charge": 0, "unpaired_electrons": 0, "electronic_temperature_K": 300.,
            "optimization_convergence": "extreme", "accuracy": .5, "energy_unit": "eV",
            "energy_description": "GFN2 total model energy including ALPB excess solvation"}
TEMPERATURES = (298.15, 360., 370., 380., 383.15, 390., 400., 410., 420.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def source_record(path):
    path = Path(path)
    return {"path": str(path.resolve()), "sha256": sha(path), "bytes": path.stat().st_size}


def single_point(atoms, folder, executable, threads, timeout):
    folder = Path(folder)
    result_path = folder / "sp_result.json"
    if result_path.exists():
        result = _json(result_path)
        if (result.get("protocol") != PROTOCOL or result.get("output_sha256") != sha(folder / "xtb.out")
                or not np.array_equal(read(folder / "input.xyz").numbers, atoms.numbers)
                or not np.allclose(read(folder / "input.xyz").positions, atoms.positions, atol=1e-8, rtol=0)):
            raise ValueError("cached single-point input, output or protocol mismatch")
        return result
    try:
        code, output, wall, argv = _execute(atoms, folder, executable, 0, 0, threads,
            ["--sp"], timeout, 300., "toluene", "gsolv")
        good = code == 0 and "normal termination of xtb" in output.lower()
        result = {"status": "converged" if good else "failed", "returncode": code,
                  "energy_eV": parse_energy(output) if good else None, "wall_seconds": wall,
                  "argv": argv, "output_sha256": sha(folder / "xtb.out"), "protocol": PROTOCOL}
    except Exception as exc:
        result = {"status": "failed", "energy_eV": None, "error": f"{type(exc).__name__}: {exc}",
                  "protocol": PROTOCOL,
                  "output_sha256": sha(folder / "xtb.out") if (folder / "xtb.out").exists() else None}
    atomic_json(result_path, result)
    return result


def checked_reference(name, smiles, out, scratch, args):
    """Reuse a source Hessian only after explicit graph, protocol and hash checks."""
    directory = REPO / "data/reference_thermochemistry" / name / "stationarity"
    record_path = directory / "result.json"
    source = _json(record_path)
    geometry = directory / "accepted_geometry.xyz"
    atoms = read(geometry)
    bonds = graph_bonds(molecular_graph(smiles))
    identity = connectivity_check(atoms, bonds)
    required = {"method": "GFN2-xTB", "solvent": "toluene", "solvation_state": "gsolv",
                "charge": 0, "unpaired_electrons": 0, "electronic_temperature_K": 300.}
    if not source.get("accepted") or not identity["retained"]:
        raise ValueError(f"{name}: source minimum or graph identity is invalid")
    if any(source.get(k) != v for k, v in required.items()) or source["accepted_geometry_sha256"] != sha(geometry):
        raise ValueError(f"{name}: source geometry hash or protocol mismatch; recalculate reference")
    attempt = source["attempts"][source["accepted_attempt"]]
    spectra = directory / f"hessian_{source['accepted_attempt']:02d}.npz"
    if any(record["source_sha256"] != sha(spectra) for record in source["thermochemistry"]):
        raise ValueError(f"{name}: reference spectrum hash mismatch")
    for key in ("native_output", "unprojected_hessian"):
        archived = attempt[key]
        if sha(directory / archived["file"]) != archived["gzip_sha256"]:
            raise ValueError(f"{name}: reference native evidence hash mismatch")
    sp = single_point(atoms, scratch / "references" / name, args.executable, args.threads, args.timeout)
    if sp["status"] != "converged":
        raise RuntimeError(f"{name}: same-protocol reference single point failed")
    record = {"species": name, "smiles": smiles, "energy_eV": sp["energy_eV"],
              "energy_hartree": sp["energy_eV"] / Hartree, "protocol": PROTOCOL,
              "source_minimum_result": source_record(record_path), "source_geometry": source_record(geometry),
              "source_spectrum": source_record(spectra), "connectivity": identity,
              "source_gradient_energy_difference_eV": sp["energy_eV"] - attempt["gradient_energy_eV"],
              "thermochemistry_reuse": "source full native Hessian, matching protocol and SHA256 verified",
              "thermochemistry": source["thermochemistry"], "new_single_point": sp}
    if abs(record["source_gradient_energy_difference_eV"]) > 1e-4:
        raise ValueError(f"{name}: fresh reference energy does not reproduce the verified source")
    atomic_json(out / "references" / f"{name}.json", record)
    write(out / "references" / f"{name}.xyz", atoms, format="xyz")
    return atoms, record


def _connected_hbond_components(atoms, fragments, bonds):
    ids = {i: f for f, fragment in enumerate(fragments) for i in fragment}
    donors = []
    for i, j in bonds:
        if atoms[i].symbol == "H" and atoms[j].symbol in ("N", "O"):
            donors.append((j, i))
        elif atoms[j].symbol == "H" and atoms[i].symbol in ("N", "O"):
            donors.append((i, j))
    edges = [{"donor": d, "hydrogen": h, "acceptor": a} for d, h in donors
             for a in range(len(atoms)) if atoms[a].symbol in ("N", "O") and ids[a] != ids[d]]
    observed = [e for e in hydrogen_bond_metrics(atoms, edges) if e["retained_after_optimization"]]
    reached = {0}
    while True:
        expanded = reached | {ids[e["donor"]] for e in observed if ids[e["acceptor"]] in reached} | {
            ids[e["acceptor"]] for e in observed if ids[e["donor"]] in reached}
        if reached == expanded:
            break
        reached = expanded
    return {"all_fragments_H_bond_connected": len(reached) == len(fragments),
            "substrate_connected_fragments": sorted(reached), "observed_hydrogen_bonds": observed}


def _csv(path, records):
    if not records:
        return
    fields = list(dict.fromkeys(k for row in records for k in row))
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def run(args):
    out, scratch = Path(args.output).resolve(), Path(args.scratch).resolve()
    out.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    summary = {"started_utc": utc_now(), "protocol": PROTOCOL, "threads": args.threads,
               "physical_workers": 1, "scope": "Open microsolvated hemiaminal minima, not TSs or barriers",
               "references": {}, "clusters": [], "selected": [], "association_thermochemistry": [],
               "n_solvent": [1, 2, 3], "orientation_seeds_per_size": args.seeds,
               "executable_source": source_record(args.executable),
               "code_sources": [source_record(__file__), source_record(REPO / "src/pincer_catmech/generators/solvation_clusters.py"),
                                source_record(REPO / "src/pincer_catmech/quantum/xtb_backend.py")]}
    substrate, subref = checked_reference("hemiaminal", HEMIAMINAL_SMILES, out, scratch, args)
    solvent, solref = checked_reference("tert_butanol", TERT_BUTANOL_SMILES, out, scratch, args)
    summary["references"] = {"hemiaminal": subref, "tert_butanol": solref}
    version = re.search(r"xtb version\s+([^\r\n]+)", (scratch / "references/hemiaminal/xtb.out").read_text(encoding="utf-8"))
    summary["program_version"] = version[1].strip() if version else "not parsed; inspect native output"
    by_name = {}
    for n in (1, 2, 3):
        for seed in range(args.seeds):
            name = f"hemiaminal_tbuoh_{n}_seed_{seed:02d}"
            folder = out / "clusters" / name
            folder.mkdir(parents=True, exist_ok=True)
            cluster = build_solvation_cluster(n, seed=seed, substrate=substrate, solvent=solvent)
            write(folder / "initial.xyz", cluster.atoms, format="xyz")
            atomic_json(folder / "mapping.json", cluster.metadata())
            native_folder = scratch / "optimizations" / name
            if (native_folder / "result.json").exists():
                native = _json(native_folder / "result.json")
                old_atoms = read(native_folder / "input.xyz")
                if (not np.array_equal(old_atoms.numbers, cluster.atoms.numbers)
                        or not np.allclose(old_atoms.positions, cluster.atoms.positions, atol=1e-8, rtol=0)):
                    raise ValueError("cached optimization geometry differs from deterministic docking")
                # The existing backend hashes decoded stdout before Windows
                # newline translation; retain a separate file-byte hash below.
                normalized_hash = hashlib.sha256((native_folder / "xtb.out").read_text(encoding="utf-8").encode()).hexdigest()
                if native.get("output_sha256") and native["output_sha256"] != normalized_hash:
                    raise ValueError("cached optimization output checksum mismatch")
                command = _json(native_folder / "command.json")
                if command.get("solvent") != "toluene" or command.get("solvation_state") != "gsolv":
                    raise ValueError("cached optimization solvent protocol mismatch")
            else:
                native = asdict(optimize_geometry(cluster.atoms, native_folder, charge=0, unpaired=0,
                    threads=args.threads, executable=args.executable, max_cycles=400,
                    timeout_seconds=args.timeout, convergence="extreme", solvent="toluene"))
            record = {"name": name, "n_solvent": n, "seed": seed, "topology": cluster.topology,
                      "native_result": native, "protocol": PROTOCOL,
                      "native_output_file": source_record(native_folder / "xtb.out"),
                      "initial_geometry": source_record(folder / "initial.xyz"),
                      "mapping": source_record(folder / "mapping.json"), "accepted_for_association": False,
                      "thermochemistry_status": "not_selected_for_hessian"}
            if native["converged"]:
                optimized = read(native_folder / "xtbopt.xyz")
                write(folder / "optimized.xyz", optimized, format="xyz")
                record["connectivity"] = connectivity_check(optimized, cluster.bonds)
                record["selected_hbond_metrics"] = hydrogen_bond_metrics(optimized, cluster.hydrogen_bonds)
                record["final_hbond_topology"] = _connected_hbond_components(optimized, cluster.fragments, cluster.bonds)
                record["accepted_for_association"] = record["connectivity"]["retained"]
                record["optimized_geometry"] = source_record(folder / "optimized.xyz")
                if record["accepted_for_association"]:
                    delta = association_energy(native["energy_eV"], subref["energy_eV"], [solref["energy_eV"]] * n)
                    record.update(association_energy_eV=delta, association_energy_kcal_mol=delta * EV_TO_KCAL_MOL,
                                  energy_hartree=native["energy_eV"] / Hartree)
                by_name[name] = (cluster, optimized)
            atomic_json(folder / "result.json", record)
            summary["clusters"].append(record)
            atomic_json(out / "summary.json", summary)
            print(json.dumps({"name": name, "status": native["status"], "association_kcal_mol": record.get("association_energy_kcal_mol"),
                              "wall_seconds": native["wall_seconds"]}), flush=True)
    for n in (1, 2, 3):
        candidates = [r for r in summary["clusters"] if r["n_solvent"] == n and r["accepted_for_association"]]
        if not candidates:
            continue
        selected = min(candidates, key=lambda r: r["native_result"]["energy_eV"])
        name = selected["name"]
        cluster, optimized = by_name[name]
        fragment_records = []
        monomers, pairs = [], {}
        failed = False
        for size in (1, 2):
            for subset in combinations(range(n + 1), size):
                indices = [i for fragment in subset for i in cluster.fragments[fragment]]
                fragment = optimized[indices]
                job = single_point(fragment, scratch / "frozen_fragments" / name / "_".join(map(str, subset)),
                                   args.executable, args.threads, args.timeout)
                fragment_records.append({"subset": list(subset), "source_atom_indices": indices, **job})
                if job["status"] != "converged":
                    failed = True
                if size == 1:
                    monomers.append(job["energy_eV"])
                else:
                    pairs[subset] = job["energy_eV"]
        selected_record = {"name": name, "n_solvent": n, "selection": "lowest energy among converged graph-retained sampled seeds",
                           "association_energy_kcal_mol": selected["association_energy_kcal_mol"],
                           "frozen_fragment_calculations": fragment_records}
        if not failed:
            full_sp = single_point(optimized, scratch / "frozen_fragments" / name / "full",
                                   args.executable, args.threads, args.timeout)
            selected_record["full_single_point"] = full_sp
            if full_sp["status"] == "converged":
                mb = many_body_nonadditivity(full_sp["energy_eV"], monomers, pairs)
                mb["deformation_energy_eV"] = sum(monomers) - subref["energy_eV"] - n * solref["energy_eV"]
                mb["beyond_pair_nonadditivity_kcal_mol"] = mb["beyond_pair_nonadditivity_eV"] * EV_TO_KCAL_MOL
                mb["interpretation"] = "Frozen-geometry beyond-pair remainder includes ALPB cavity and reaction-field nonadditivity"
                selected_record["many_body"] = mb
        if args.hessians and selected["final_hbond_topology"]["all_fragments_H_bond_connected"]:
            thermo_dir = out / "clusters" / name / "stationarity"
            if (thermo_dir / "result.json").exists():
                thermo = _json(thermo_dir / "result.json")
                original_map = thermo.get("source_metadata", {}).get("atom_mapping")
                if original_map and original_map["sha256"] != selected["mapping"]["sha256"]:
                    # Closure diagnostics were added after the physical run.
                    # Reconstruct and preserve the byte-identical original map,
                    # validated against the original recorded SHA256.
                    prior_mapping = thermo_dir.parent / "mapping_before_closure_probes.json"
                    if not prior_mapping.exists():
                        prior = _json(thermo_dir.parent / "mapping.json")
                        prior.pop("initial_closure_probes", None)
                        atomic_json(prior_mapping, prior)
                    if sha(prior_mapping) != original_map["sha256"]:
                        raise ValueError("cannot reconstruct provenance-tracked original atom mapping")
                    thermo["source_metadata"]["atom_mapping"] = source_record(prior_mapping)
                    thermo["mapping_source_relocation_audit"] = {
                        "operation": "preserve original source mapping after adding closure diagnostics",
                        "original_source_sha256": original_map["sha256"], "content_changed": False}
                    atomic_json(thermo_dir / "result.json", thermo)
                if thermo["accepted"] and tuple(thermo["temperatures_K"]) != TEMPERATURES:
                    # Extend temperature evaluation of the same certified Hessian.
                    # Preserve the first physical-run record before metadata updates.
                    history = thermo_dir / "metadata_history"
                    history.mkdir(exist_ok=True)
                    old_record = history / "result_before_temperature_extension.json"
                    if not old_record.exists():
                        shutil.copy2(thermo_dir / "result.json", old_record)
                    accepted_atoms = read(thermo_dir / "accepted_geometry.xyz")
                    attempt = thermo["attempts"][thermo["accepted_attempt"]]
                    spectrum = thermo_dir / f"hessian_{thermo['accepted_attempt']:02d}.npz"
                    if (sha(thermo_dir / "accepted_geometry.xyz") != thermo["accepted_geometry_sha256"] or
                            any(t["source_sha256"] != sha(spectrum) for t in thermo["thermochemistry"])):
                        raise ValueError("cached thermochemical geometry or Hessian checksum mismatch")
                    thermo["thermochemistry"] = [asdict(r) for r in thermochemistry_from_modes(
                        accepted_atoms, attempt["energy_eV"], attempt["frequencies_cm1"],
                        temperatures=TEMPERATURES, unpaired=0, symmetry_number=1, solvent="toluene",
                        solvation_state="gsolv", source=str(spectrum), source_sha256=sha(spectrum))]
                    thermo["temperatures_K"] = list(TEMPERATURES)
                    thermo["temperature_extension_audit"] = {"utc": utc_now(), "prior_result": source_record(old_record),
                        "operation": "re-evaluate identical certified native Hessian at additional temperatures",
                        "new_native_calculations": 0, "energies_and_geometries_modified": False}
                    atomic_json(thermo_dir / "result.json", thermo)
            else:
                thermo = characterize_stationary_point(optimized, thermo_dir, scratch / "stationarity" / name,
                    charge=0, unpaired=0, threads=args.threads, executable=args.executable,
                    max_repairs=1, temperatures=TEMPERATURES, solvent="toluene",
                    deadline_epoch=time.time() + args.hessian_timeout,
                    source_metadata={"name": name, "n_solvent": n, "optimization": selected["optimized_geometry"],
                                     "atom_mapping": selected["mapping"], "state_model": "isolated microsolvated molecular cluster"})
            selected["thermochemistry_status"] = thermo["status"]
            selected_record["thermochemistry_status"] = thermo["status"]
            selected_record["stationarity_result"] = source_record(thermo_dir / "result.json")
            if thermo["accepted"]:
                accepted = read(thermo_dir / "accepted_geometry.xyz")
                graph_ok = connectivity_check(accepted, cluster.bonds)["retained"]
                selected_record["thermochemical_geometry_connectivity_retained"] = graph_ok
                selected_record["thermochemical_hbond_topology"] = _connected_hbond_components(accepted, cluster.fragments, cluster.bonds)
                if graph_ok and selected_record["thermochemical_hbond_topology"]["all_fragments_H_bond_connected"]:
                    for row in thermo["thermochemistry"]:
                        temp = row["temperature"]
                        sr = next(x for x in subref["thermochemistry"] if x["temperature"] == temp)
                        ar = next(x for x in solref["thermochemistry"] if x["temperature"] == temp)
                        dg = row["G_298_qRRHO_sol"] - sr["G_298_qRRHO_sol"] - n * ar["G_298_qRRHO_sol"]
                        summary["association_thermochemistry"].append({"name": name, "n_solvent": n,
                            "temperature_K": temp, "delta_G_qRRHO_1M_kcal_mol": dg,
                            "standard_state": "all species 1 M; ALPB gsolv plus single 1 atm to 1 M correction",
                            "status": "verified_cluster_and_reference_minima",
                            "interpretation": "sampled cluster association; not dehydration barrier or bulk speciation"})
        else:
            selected_record["thermochemistry_status"] = "disabled" if not args.hessians else "fragment_hbond_disconnection"
        summary["selected"].append(selected_record)
        atomic_json(out / "clusters" / name / "result.json", selected)
        atomic_json(out / "summary.json", summary)
        print(json.dumps({"selected": name, "nonadditivity_kcal_mol": selected_record.get("many_body", {}).get("beyond_pair_nonadditivity_kcal_mol"),
                          "thermochemistry_status": selected_record["thermochemistry_status"]}), flush=True)
    rows = [{"name": r["name"], "n_solvent": r["n_solvent"], "seed": r["seed"],
             "status": r["native_result"]["status"], "connectivity_retained": r.get("connectivity", {}).get("retained"),
             "all_fragments_H_bond_connected": r.get("final_hbond_topology", {}).get("all_fragments_H_bond_connected"),
             "energy_eV": r["native_result"]["energy_eV"], "association_energy_kcal_mol": r.get("association_energy_kcal_mol"),
             "max_force_eV_A": r["native_result"]["max_force_eV_A"], "thermochemistry_status": r["thermochemistry_status"]}
            for r in summary["clusters"]]
    _csv(out / "association_energies.csv", rows)
    _csv(out / "association_thermochemistry.csv", summary["association_thermochemistry"])
    _csv(out / "many_body_nonadditivity.csv", [{"name": r["name"], "n_solvent": r["n_solvent"],
          **{k: v for k, v in r.get("many_body", {}).items() if not isinstance(v, dict)}} for r in summary["selected"]])
    manifest = []
    with zipfile.ZipFile(out / "native_evidence.zip", "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(scratch.rglob("*")):
            if path.is_file() and path.name != ".worker.lock":
                relative = path.relative_to(scratch).as_posix()
                archive.write(path, relative)
                manifest.append({"path": relative, "sha256": sha(path), "bytes": path.stat().st_size})
    atomic_json(out / "native_manifest.json", manifest)
    summary.update(finished_utc=utc_now(), native_archive=source_record(out / "native_evidence.zip"),
                   native_manifest=source_record(out / "native_manifest.json"),
                   converged_optimizations=sum(r["native_result"]["converged"] for r in summary["clusters"]),
                   limitation="Bounded conformer/motif sampling. No tBuOK ion pair, catalyst, TS/IRC or solution population determination.")
    atomic_json(out / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPO / "data/phase3/solvation")
    parser.add_argument("--scratch", type=Path, default=REPO.parent / "phase3/solvation")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--executable")
    parser.add_argument("--timeout", type=float, default=300.)
    parser.add_argument("--hessian-timeout", type=float, default=420.)
    parser.add_argument("--hessians", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    if not 1 <= args.seeds <= 8 or not 1 <= args.threads <= 2 or min(args.timeout, args.hessian_timeout) <= 0:
        parser.error("bounded run requires 1..8 seeds, 1..2 threads and positive time limits")
    args.executable = executable_path(args.executable)
    args.scratch.mkdir(parents=True, exist_ok=True)
    lock = args.scratch / ".worker.lock"
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(str(os.getpid()))
        result = run(args)
        print(json.dumps({"finished_utc": result["finished_utc"], "converged": result["converged_optimizations"],
                          "thermochemistry_rows": len(result["association_thermochemistry"])}), flush=True)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
