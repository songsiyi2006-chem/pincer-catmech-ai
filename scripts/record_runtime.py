"""Record the actual scientific runtime without importing large model libraries."""
from importlib import metadata
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import psutil

REPO = Path(__file__).resolve().parents[1]
PACKAGES = ("numpy", "scipy", "ase", "rdkit", "scikit-learn", "psutil", "pytest",
            "pydantic", "reportlab", "matplotlib", "pillow", "jieba")


def main():
    executable = Path(os.environ["PINCER_XTB"]).resolve()
    hasher = hashlib.sha256()
    with executable.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    versions = {}
    for name in PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    command = subprocess.run([str(executable), "--version"], capture_output=True, text=True,
                             encoding="utf-8", errors="replace", timeout=30)
    if command.returncode:
        raise RuntimeError(command.stderr)
    value = {"python": sys.version, "python_executable": sys.executable,
        "platform": platform.platform(), "processor": platform.processor(),
        "logical_CPUs": psutil.cpu_count(), "physical_cores": psutil.cpu_count(logical=False),
        "physical_memory_bytes": psutil.virtual_memory().total,
        "computational_device": "CPU; native xTB OMP resources are recorded per job",
        "packages": versions, "native_xtb_executable": str(executable),
        "native_xtb_sha256": hasher.hexdigest(), "native_xtb_version_stdout": command.stdout,
        "native_xtb_version_stderr": command.stderr,
        "existing_environment_preserved": "Campaign overlay venv inherits phase2ff; only missing test/report helpers installed in the overlay",
        "hosted_language_model_CPU": "Not observable; telemetry measures local scientific subprocesses only"}
    folder = REPO / "data/verification"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "RUNTIME.json").write_text(json.dumps(value, indent=2), encoding="utf-8")
    (folder / "requirements-observed.txt").write_text(
        "# Observed Python distributions in this campaign; not a cross-platform lockfile.\n"
        "# Native xTB and its binary hash are recorded separately in RUNTIME.json.\n" +
        "".join(f"{name}=={version}\n" for name, version in versions.items() if version), encoding="utf-8")
    print(json.dumps({"packages": versions, "logical_CPUs": value["logical_CPUs"],
                      "physical_memory_bytes": value["physical_memory_bytes"], "xtb_sha256": hasher.hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
