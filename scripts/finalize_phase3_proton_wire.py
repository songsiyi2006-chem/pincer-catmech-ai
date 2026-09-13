"""Archive, verify and tabulate completed proton-wire runs without QC calls.

New archives preserve every raw scratch file, split by run and compressed size.
An already finalized directory is verified read-only; its receipt is unchanged.
"""
from pathlib import Path, PurePosixPath
import argparse
import csv
import hashlib
import json
import sys
import uuid
import zipfile
import zlib

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from pincer_catmech.quantum.xtb_backend import atomic_json, utc_now

MAX_ARCHIVE_BYTES = 95 * 1024 * 1024


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def separate_paths(output, scratch):
    output, scratch = Path(output).resolve(), Path(scratch).resolve()
    if output == scratch or output.is_relative_to(scratch) or scratch.is_relative_to(output):
        raise ValueError("Output and scratch must be separate, non-nested directories")
    return output, scratch


def checked_member(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError("Unsafe archive member path")
    return name


def recorded_name(value):
    """Read portable basenames from native Windows or POSIX provenance."""
    return checked_member(PurePosixPath(str(value).replace("\\", "/")).name)


def archive_native_tree(scratch, output, *, max_archive_bytes=MAX_ARCHIVE_BYTES):
    """Write immutable parts; sizing includes ZIP central-directory overhead."""
    output, scratch = separate_paths(output, scratch)
    if not 512 <= max_archive_bytes <= MAX_ARCHIVE_BYTES:
        raise ValueError("Archive limit must be between 512 bytes and 95 MiB")
    if not scratch.is_dir():
        raise FileNotFoundError(scratch)
    if (output / "native_manifest.json").exists() or list(output.glob("native_evidence_*.zip")):
        raise FileExistsError("Existing native archives are immutable")
    output.mkdir(parents=True, exist_ok=True)
    staging = output / (".native_staging_" + uuid.uuid4().hex)
    staging.mkdir()
    members, part_paths, writer = [], [], None
    planned_bytes, member_count, active_group = 98, 0, None
    try:
        for path in sorted(scratch.rglob("*")):
            if not path.is_file():
                continue
            if path.is_symlink() or not path.resolve().is_relative_to(scratch):
                raise ValueError("Native evidence links may not escape scratch")
            member = checked_member(path.relative_to(scratch).as_posix())
            payload = path.read_bytes()
            compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
            compressed_size = len(compressor.compress(payload) + compressor.flush())
            # Local and central headers, UTF-8 member names, plus ZIP64 reserve.
            addition = compressed_size + 76 + 2 * len(member.encode("utf-8")) + 40
            if addition + 98 >= max_archive_bytes:
                raise ValueError(f"Single compressed member exceeds archive limit: {member}")
            group = PurePosixPath(member).parts[0]
            if writer is None or group != active_group or planned_bytes + addition >= max_archive_bytes or member_count >= 65000:
                if writer is not None:
                    writer.close()
                part = staging / f"native_evidence_{len(part_paths)+1:03d}.zip"
                part_paths.append(part)
                writer = zipfile.ZipFile(part, "w", zipfile.ZIP_DEFLATED, compresslevel=6)
                planned_bytes, member_count, active_group = 98, 0, group
            writer.writestr(member, payload)
            members.append({"archive": part.name, "member": member,
                "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)})
            planned_bytes += addition
            member_count += 1
    finally:
        if writer is not None:
            writer.close()
    archives = []
    for part in part_paths:
        if part.stat().st_size >= max_archive_bytes:
            raise ValueError("Final ZIP size exceeds the explicit maximum; raw files remain intact")
        archives.append({"name": part.name, "sha256": sha256(part), "bytes": part.stat().st_size,
            "member_count": sum(row["archive"] == part.name for row in members)})
    # Publish only after every part meets the size bound; never replace a part.
    for part in part_paths:
        destination = output / part.name
        if destination.exists():
            raise FileExistsError(destination)
        part.rename(destination)
    staging.rmdir()
    manifest = {"schema": "native_evidence_multiarchive_v1", "archives": archives, "members": members,
        "member_count": len(members), "unique_current_member_count": len(members),
        "raw_bytes": sum(row["bytes"] for row in members), "archive_max_bytes": max_archive_bytes,
        "all_current_members_covered": True, "scratch_source": str(scratch),
        "scope": "Every raw scratch file, including electronic restart/cache files; one archive-qualified copy per file"}
    atomic_json(output / "native_manifest.json", manifest)
    return manifest


def verify_archives(output):
    """Verify archive-qualified versions, including the historical two-part run."""
    output = Path(output).resolve()
    manifest = json.loads((output / "native_manifest.json").read_text(encoding="utf-8"))
    archives = manifest.get("archives")
    if archives is None:
        archives = [{"name": "native_evidence.zip", **manifest["archive"]}]
    seen, verified = set(), 0
    for item in archives:
        name = checked_member(item["name"])
        if PurePosixPath(name).name != name or name in seen:
            raise ValueError("Archive names must be unique local basenames")
        seen.add(name)
        path = output / name
        if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError("Archive SHA256 or byte count mismatch: " + name)
        if path.stat().st_size >= MAX_ARCHIVE_BYTES:
            raise ValueError("Archive is not below 95 MiB")
        expected = [row for row in manifest["members"] if row.get("archive", name) == name]
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            if len(names) != len(expected) or len(set(names)) != len(names) or set(names) != {row["member"] for row in expected}:
                raise ValueError("Archive membership mismatch")
            for row in expected:
                payload = archive.read(checked_member(row["member"]))
                if len(payload) != row["bytes"] or hashlib.sha256(payload).hexdigest() != row["sha256"]:
                    raise ValueError("Member SHA256 or byte count mismatch")
                verified += 1
    if verified != len(manifest["members"]):
        raise ValueError("Manifest contains members without a matching archive")
    return {"archive_members_verified": verified, "archives_verified": archives}


def tabulate_paths(output, summary):
    """Inspect recorded coordinates/energies only; never attach a calculator."""
    import numpy as np
    from ase.io import read
    from pincer_catmech.kinetics.condensation_search import identity

    bands, endpoints, thermochemistry, attempts = [], [], [], []
    for attempt in summary.get("attempts", []):
        # Resolve by basename so a copied output tree remains auditable.
        path = output / recorded_name(attempt["directory"])
        result_path = path / "result.json"
        if sha256(result_path) != attempt["result"]["sha256"]:
            raise ValueError("Recorded path result SHA256 differs")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        protocol = json.loads((path / "protocol.json").read_text(encoding="utf-8"))
        if protocol["threads"] != 1 or protocol["charge"] != 0 or protocol["solvent"] != "toluene" or not protocol["full_endpoint_hessians_required"]:
            raise ValueError("Incompatible neutral proton-wire protocol")
        mapping = json.loads((output / (path.name + "_assembly") / "mapping_and_provenance.json").read_text(encoding="utf-8"))
        indices = mapping["mapped_indices"]
        if "restart" in result:
            restart = result["restart"]
            portable = PurePosixPath(str(restart["source"]).replace("\\", "/"))
            restart_path = output / checked_member(portable.parent.name) / checked_member(portable.name)
            if not restart["all_forces_recomputed"] or restart["new_interpolation"] or sha256(restart_path) != restart["sha256"]:
                raise ValueError("Restart provenance or recomputation gate failed")
            original = read(restart_path, index=":")
            initial_path = path / "initial_band.xyz"
            if initial_path.exists():
                initial = read(initial_path, index=":")
                if len(original) != 9 or len(initial) != 9:
                    raise ValueError("Nine restart images are required")
                for before, after in zip(original, initial):
                    np.testing.assert_array_equal(before.numbers, after.numbers)
                    np.testing.assert_allclose(before.positions, after.positions, rtol=0, atol=1e-8)
        logs = (path / "neb.log").read_text().splitlines() if (path / "neb.log").exists() else []
        fire = [line for line in logs if line.startswith("FIRE:")]
        attempts.append({"attempt": path.name, "status": result["status"], "accepted_TS": attempt["accepted_TS"],
            "last_NEB_step": int(fire[-1].split()[1]) if fire else attempt["last_reported_NEB_step"],
            "last_NEB_fmax_eV_A": float(fire[-1].split()[-1]) if fire else None,
            "NEB_fmax_target_eV_A": protocol["neb_fmax_eV_A"],
            "last_band_max_above_reactant_eV": attempt["last_sampled_band_potential_max_above_reactant_eV"],
            "is_activation_free_energy": False,
            "original_product_preparation_graph_retained": mapping["product_seed_identity"]["retained"]})
        for endpoint in result.get("endpoints", []):
            index = endpoint["index"]
            frequencies = endpoint.get("frequencies_cm1", [])
            if endpoint.get("stationarity_status") == "verified_minimum":
                geometry = read(path / f"endpoint_{index}.xyz")
                graph = identity(geometry, protocol["endpoint_bonds"][index], cn=protocol["C_N_indices"],
                    cn_interval=(1.35, 1.70) if index == 0 else (1.15, 1.42))
                if not graph["retained"] or len(frequencies) != 3*len(geometry)-6 or min(frequencies) <= 0:
                    raise ValueError("Endpoint graph or full positive-frequency count failed")
                if sha256(path / f"endpoint_{index}_hessian.npz") != endpoint["hessian_sha256"]:
                    raise ValueError("Endpoint Hessian SHA256 differs")
            endpoints.append({"attempt": path.name, "endpoint_index": index, "stationarity_status": endpoint.get("stationarity_status"),
                "graph_retained": endpoint["identity"]["retained"], "C_N_distance_A": endpoint["identity"]["C_N_distance_angstrom"],
                "minimum_frequency_cm1": min(frequencies) if frequencies else None,
                "imaginary_modes": sum(f < 0 for f in frequencies), "electronic_energy_eV": endpoint["energy_eV"]})
            for row in endpoint.get("thermochemistry", []):
                thermochemistry.append({"attempt": path.name, "endpoint_index": index, "temperature_K": row["temperature"],
                    "E_elec_kcal_mol": row["E_elec"], "ZPVE_kcal_mol": row["ZPVE"], "G_qRRHO_1M_kcal_mol": row["G_298_qRRHO_sol"],
                    "S_qRRHO_cal_mol_K": row["S_qRRHO"], "cutoff_cm1": row["cutoff_cm1"], "is_activation_free_energy": False})
        band_path = path / "final_band.xyz"
        if not band_path.exists():
            band_path = path / "current_band.xyz"
        if band_path.exists():
            frames = read(band_path, index=":")
            snapshots = [row for row in result["stages"] if row["stage"] == "neb"]
            energies = result.get("image_energies_eV") or (snapshots[-1]["image_energies_eV"] if snapshots else [])
            if len(frames) != 9 or len(energies) != 9:
                raise ValueError("A recorded band must have nine geometries and energies")
            for i, (atoms, energy) in enumerate(zip(frames, energies)):
                np.testing.assert_array_equal(atoms.numbers, frames[0].numbers)
                if len(atoms) != 43:
                    raise ValueError("Expected a 43-atom mapped neutral cluster")
                row = {"attempt": path.name, "image": i, "energy_eV": energy, "relative_potential_eV": energy-energies[0],
                    "relative_potential_kcal_mol": (energy-energies[0])*23.060547830619,
                    "band_converged": attempt["neb_converged"], "is_activation_free_energy": False}
                for label, pair in {"C_O_A": ("C", "leaving_O"), "C_N_A": ("C", "N"), "N_H_A": ("N", "N_H"),
                    "tBuOH_O_original_N_H_A": ("shuttle_O", "N_H"), "tBuOH_O_original_alcohol_H_A": ("shuttle_O", "shuttle_H"),
                    "water_O_original_alcohol_H_A": ("leaving_O", "shuttle_H")}.items():
                    row[label] = atoms.get_distance(indices[pair[0]], indices[pair[1]])
                bands.append(row)
    for filename, rows in (("sampled_band_profiles.csv", bands), ("endpoint_stationarity.csv", endpoints),
                           ("endpoint_thermochemistry.csv", thermochemistry), ("attempts.csv", attempts)):
        if rows:
            with (output / filename).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    return {"completed_path_runs": len(attempts), "distinct_initial_seed_attempts": sum("continuation" not in row["attempt"] for row in attempts),
        "continuation_count": sum("continuation" in row["attempt"] for row in attempts),
        "verified_endpoint_minimum_certifications": sum(row["stationarity_status"] == "verified_minimum" for row in endpoints),
        "initial_seed_endpoint_certifications": sum(row["stationarity_status"] == "verified_minimum" and "continuation" not in row["attempt"] for row in endpoints),
        "actual_band_image_rows": len(bands), "endpoint_thermochemistry_rows": len(thermochemistry), "attempts": attempts}


def finalize_evidence(output, scratch, *, max_archive_bytes=MAX_ARCHIVE_BYTES):
    output, scratch = separate_paths(output, scratch)
    summary_path = output / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary["status"] == "running":
        raise ValueError("Wait for physical work to finish before finalization")
    if (output / "verification.json").exists():
        receipt = json.loads((output / "verification.json").read_text(encoding="utf-8"))
        if receipt.get("status") != "passed":
            raise ValueError("Existing non-passed receipt requires explicit inspection")
        return {"already_finalized": True, "existing_receipt_unchanged": True, **verify_archives(output)}
    if not (output / "native_manifest.json").exists():
        archive_native_tree(scratch, output, max_archive_bytes=max_archive_bytes)
    audit = verify_archives(output)
    original = output / "finalization_input_summary.json"
    if not original.exists():
        original.write_bytes(summary_path.read_bytes())
    tables = tabulate_paths(output, summary)
    sources = summary.get("executed_driver_sources", [])
    for record in sources:
        snapshot = output / checked_member(record["relative_path"])
        if sha256(snapshot) != record["sha256"]:
            raise ValueError("Executed driver snapshot SHA256 differs")
    source_snapshot = output / "evidence_finalizer_source.py"
    if not source_snapshot.exists():
        source_snapshot.write_bytes(Path(__file__).read_bytes())
    receipt = {"verified_utc": utc_now(), "status": "passed", **audit, **tables,
        "accepted_TS_count": summary["accepted_TS_count"], "activation_free_energy_barrier_count": summary["activation_free_energy_barrier_count"],
        "baseline_temperature_K": summary.get("thermodynamic_baseline_temperature_K"),
        "executed_driver_sources": sources, "executed_driver_source_status": "snapshots_verified" if sources else "not_recorded_by_this_run",
        "finalizer_source_sha256": sha256(source_snapshot), "new_quantum_calculations_performed": False,
        "limits": "Band maxima are sampled potential diagnostics; TS and free-energy-barrier acceptance comes only from the recorded physical gates."}
    manifest = json.loads((output / "native_manifest.json").read_text(encoding="utf-8"))
    summary["native_evidence"] = {"schema": manifest["schema"], "archives": manifest["archives"],
        "member_count": manifest["member_count"], "unique_current_member_count": manifest["unique_current_member_count"],
        "all_current_members_covered": manifest["all_current_members_covered"]}
    atomic_json(summary_path, summary)
    atomic_json(output / "verification.json", receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--verify-only", "--verify-evidence", dest="verify_only", action="store_true")
    args = parser.parse_args()
    if args.verify_only:
        result = verify_archives(args.output)
    else:
        if args.scratch is None:
            parser.error("--scratch is required for finalization")
        result = finalize_evidence(args.output, args.scratch)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
