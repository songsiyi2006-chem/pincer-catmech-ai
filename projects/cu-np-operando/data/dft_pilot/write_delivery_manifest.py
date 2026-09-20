"""List only deliverable files; large rebuildable scratch remains local."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
files = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file():
        continue
    relative = path.relative_to(ROOT)
    if "scratch" in relative.parts or "__pycache__" in relative.parts or path.name.startswith("psi.") or path.suffix == ".pyc" or path.name == "delivery_files.json":
        continue
    files.append({"path": relative.as_posix(), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
data = {"files": files, "file_count": len(files), "total_bytes": sum(r["bytes"] for r in files), "policy": "Copy these files retaining relative paths, plus this manifest; do not copy scratch/, __pycache__/, *.pyc or root psi.*"}
(ROOT / "delivery_files.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: data[k] for k in ["file_count", "total_bytes"]}))
