"""Archive native scientific evidence without installed software or restart caches.

Archive members are relative to --work-root. Restore them there to resolve the
historical run layout; derived datasets use portable repository-relative paths.
Every included byte is hashed, archived, and checked by reading the ZIP back.
Run only after scientific writers finish. --inventory performs no writes.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import zipfile

REPO = Path(__file__).resolve().parents[1]
FAMILIES = ("campaign-runs", "campaign-runs-alpb-toluene", "neb-scratch",
            "reference-scratch", "targeted-rescue-scratch", "cobalt-carbonyl-scratch",
            "ionic-reference-scratch", "hessian-runs", "hessian-runs-identity-replacements",
            "condensation-scratch", "co-neb-retry-scratch", "quantum-pilot")
NAMES = {"gradient", "hessian", "wbo", "charges", "energy", "vibspectrum", "xtbhess.coord"}
SUFFIXES = {".xyz", ".out", ".log", ".json", ".jsonl", ".inp", ".in", ".txt", ".engrad"}


def selected(path):
    return path.name in NAMES or path.suffix.lower() in SUFFIXES


def enumerate_files(work_root, families=FAMILIES):
    for family in families:
        # Directory walking avoids materializing every binary cache path before
        # selecting the scientific records in these large native work trees.
        for folder, directories, names in os.walk(work_root / family):
            directories.sort()
            for name in sorted(names):
                path = Path(folder) / name
                if selected(path) and path.is_file():
                    yield family, path


def package(work_root, output, *, inventory=False, target_bytes=128 * 1024**2, families=FAMILIES):
    if not families or any(family not in FAMILIES for family in families):
        raise ValueError("Select known scientific scratch families")
    files = list(enumerate_files(work_root, families))
    totals = Counter()
    for family, path in files:
        totals[family] += path.stat().st_size
    if inventory:
        return {"files": len(files), "uncompressed_bytes_by_family": dict(totals)}
    output.mkdir(parents=True, exist_ok=True)
    if list(output.glob("native-*.zip")):
        raise FileExistsError("Native archives exist; use a fresh output directory")
    entries, archives = [], []
    for family in families:
        groups, current, size = [], [], 0
        for current_family, path in files:
            if current_family != family:
                continue
            length = path.stat().st_size
            if current and size + length > target_bytes:
                groups.append(current)
                current, size = [], 0
            current.append(path)
            size += length
        if current:
            groups.append(current)
        for number, group in enumerate(groups, 1):
            archive = output / f"native-{family}-{number:03d}.zip"
            members = []
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
                for path in group:
                    content = path.read_bytes()
                    member = path.relative_to(work_root).as_posix()
                    bundle.writestr(member, content)
                    entry = {"archive": archive.name, "member": member, "bytes": len(content),
                             "sha256": hashlib.sha256(content).hexdigest()}
                    members.append(entry)
                    entries.append(entry)
            with zipfile.ZipFile(archive) as bundle:
                if bundle.testzip() is not None:
                    raise ValueError(f"CRC failure: {archive}")
                for entry in members:
                    if hashlib.sha256(bundle.read(entry["member"])).hexdigest() != entry["sha256"]:
                        raise ValueError(f"Archive hash failure: {entry['member']}")
            if archive.stat().st_size >= 95 * 1024**2:
                raise ValueError("Archive exceeds publication size bound")
            archives.append({"file": archive.name, "bytes": archive.stat().st_size,
                             "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                             "verified_members": len(members)})
            print(json.dumps({"verified_archive": archive.name, "members": len(members),
                              "bytes": archive.stat().st_size}), flush=True)
    manifest = {"schema_version": 1, "archive_format": "ZIP/deflate; work-root-relative members",
        "original_work_root": str(work_root), "families": list(families),
        "selection": {"exact_names": sorted(NAMES), "suffixes": sorted(SUFFIXES)},
        "omitted": "Installed runtimes, executables, binary restart/cache files, xtbtopo.mol topology caches, .xtboptok sentinels and all names/extensions outside the explicit selection. Restart provenance and hashes remain in launch JSON; energies/forces/geometries/native outputs are included. Molecular identity is audited with saved XYZ, expected graph metadata, native WBO and full output, not the omitted topology cache.",
        "integrity": "All archived members read back; CRC and full SHA256 verified",
        "archives": archives, "files": entries, "file_count": len(entries),
        "uncompressed_bytes_by_family": dict(totals)}
    (output / "NATIVE_EVIDENCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {"file_count": len(entries), "archive_count": len(archives),
            "compressed_bytes": sum(a["bytes"] for a in archives)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, default=REPO.parent)
    parser.add_argument("--output", type=Path, default=REPO / "data/native_evidence")
    parser.add_argument("--inventory", action="store_true")
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=FAMILIES)
    args = parser.parse_args()
    print(json.dumps(package(args.work_root.resolve(), args.output.resolve(), inventory=args.inventory, families=args.families), indent=2))
