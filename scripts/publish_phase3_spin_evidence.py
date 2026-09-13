"""Index published spin evidence and locally retained large Psi4 intermediates.

Never deletes, moves, compresses or edits native files. Codes 97 (DF-SCF B matrix)
and 64 (DIIS storage) are excluded from Git only, by explicit .gitignore rules.
All input/output/gradient/result files and other numbered files remain published.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]
RECEIPT = "evidence_publication.json"


def cache_kind(path):
    match = re.fullmatch(r"psi\.\d+\.(97|64)", Path(path).name)
    return {"97": "PSIF_DFSCF_BJ", "64": "PSIF_LIBDIIS"}.get(match[1]) if match else None


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def publish(directory):
    directory = Path(directory).resolve()
    state = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    if not state.get("run_complete") or state.get("target_states_recorded") != 54:
        raise RuntimeError("Finish the 54-slot spin matrix before indexing native evidence")
    retained, cached = [], []
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path == directory / RECEIPT:
            continue
        if path.suffix == ".tmp":
            raise RuntimeError("Unfinished atomic write found; wait for native runners")
        before = path.stat()
        entry = {"path": path.relative_to(directory).as_posix(), "bytes": before.st_size, "sha256": sha(path)}
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError("Native file changed during indexing")
        kind = cache_kind(path)
        if kind:
            entry.update(kind=kind, retained_locally=True, local_path=str(path), published_in_git=False)
            cached.append(entry)
        else:
            if before.st_size >= 100 * 1024**2:
                raise RuntimeError(f"Scientific artifact needs a lossless split archive before Git publication: {path}")
            retained.append(entry)
    result = {"schema": "phase3_spin_publication_v1", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "published_files": retained, "local_only_caches": cached,
        "published_bytes": sum(x["bytes"] for x in retained), "local_cache_bytes": sum(x["bytes"] for x in cached),
        "source": "https://psi4.github.io/psi4docs/master/autodoc_psifiles.html",
        "policy": "Only numbered DF-SCF B-matrix (97) and DIIS (64) intermediates are omitted from Git; all remain unchanged on the originating workstation. This is not a full scratch backup.",
        "scientific_evidence": "Inputs, commands, stdout, Psi4 output, analytic gradients, SCF/state responses, failures, geometries and source hashes are published."}
    (directory / RECEIPT).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def verify(directory, *, local_caches=False):
    directory = Path(directory).resolve()
    result = json.loads((directory / RECEIPT).read_text(encoding="utf-8"))
    if result.get("schema") != "phase3_spin_publication_v1":
        raise RuntimeError("Unknown spin publication manifest")
    all_entries = result["published_files"] + result["local_only_caches"]
    if len(all_entries) != len({x["path"] for x in all_entries}):
        raise RuntimeError("Duplicate publication path")
    for category in ("published_files", "local_only_caches"):
        for entry in result[category]:
            path = (directory / entry["path"]).resolve()
            if not path.is_relative_to(directory):
                raise RuntimeError("Publication path outside spin directory")
            if category == "local_only_caches":
                if cache_kind(path) != entry["kind"]:
                    raise RuntimeError("Scientific evidence cannot be relabeled as omitted cache")
                if not local_caches:
                    continue
            if sha(path) != entry["sha256"] or path.stat().st_size != entry["bytes"]:
                raise RuntimeError(f"Spin evidence SHA/size mismatch: {entry['path']}")
    observed = {p.relative_to(directory).as_posix() for p in directory.rglob("*")
                if p.is_file() and p.name != RECEIPT and not cache_kind(p)}
    if observed != {x["path"] for x in result["published_files"]}:
        raise RuntimeError("Published spin file membership differs")
    if local_caches:
        observed_caches = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() and cache_kind(p)}
        if observed_caches != {x["path"] for x in result["local_only_caches"]}:
            raise RuntimeError("Local cache inventory differs")
    return {"published_files_verified": len(result["published_files"]),
        "cache_files_retained_locally": len(result["local_only_caches"]),
        "local_cache_bytes": result["local_cache_bytes"], "local_caches_verified": local_caches,
        "publication_manifest_sha256": sha(directory / RECEIPT), "verified": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=REPO / "data/phase3/spin")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--verify-local-caches", action="store_true")
    args = parser.parse_args()
    if not args.verify:
        publish(args.directory)
    print(json.dumps(verify(args.directory, local_caches=args.verify_local_caches)), flush=True)


if __name__ == "__main__":
    main()
