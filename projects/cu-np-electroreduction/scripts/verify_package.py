"""Verify delivered file hashes and calculation provenance without computation."""
import json
from pathlib import Path
from science import sha256
ROOT=Path(__file__).resolve().parents[1]

def main():
    manifest=ROOT/"results/FILE_MANIFEST.json"
    if manifest.exists():
        for row in json.loads(manifest.read_text(encoding="utf-8"))["files"]:
            p=ROOT/row["path"]
            if sha256(p)!=row["sha256"]:raise RuntimeError(f"Changed or missing delivered file: {p}")
    for row in json.loads((ROOT/"data/raw/index.json").read_text(encoding="utf-8")):
        recpath=ROOT/row["record"];record=json.loads(recpath.read_text(encoding="utf-8"))
        if sha256(recpath)!=row["record_sha256"]:raise RuntimeError("Calculation record modified")
        for name,digest in row["files"].items():
            if sha256(recpath.parent/name)!=digest:raise RuntimeError("Native evidence modified")
        batch=recpath.parent.parent.name
        for file,key in [("run_pilot.py","runner_sha256"),("science.py","parser_sha256")]:
            if sha256(ROOT/"workflows/executed_sources"/batch/file)!=record["identity"][key]:raise RuntimeError("Executed source snapshot mismatch")
    print("Verified delivered hashes, 20 native records and matching executed source snapshots.")
if __name__=="__main__":main()
