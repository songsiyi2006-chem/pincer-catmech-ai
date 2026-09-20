"""Run one Psi4 child with wall-time and sampled process-tree memory limits.

The runner needs psutil; the child interpreter needs Psi4. No scheduler, cloud
service or software installation is invoked. Existing job output is immutable.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import psutil


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi4-python", required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--method", choices=["pbe0", "b3lyp"], required=True)
    ap.add_argument("--charge", type=int, choices=[0, -1], required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = args.output.resolve()
    if root.exists():
        raise SystemExit("Refusing to overwrite an existing calculation directory")
    available = psutil.virtual_memory().available
    if available < 2 * 1024**3:
        raise SystemExit(f"Insufficient free RAM for <=1 GiB child limit: {available} bytes available")
    root.mkdir(parents=True)
    shutil.copy2(args.geometry, root / "input.xyz")
    script = Path(__file__).with_name("psi4_single_point.py")
    shutil.copy2(script, root / "executed_psi4_single_point.py")
    shutil.copy2(__file__, root / "executed_bounded_runner.py")
    cfg = {
        "method": args.method, "basis": "def2-svpd", "charge": args.charge,
        "multiplicity": 1 if args.charge == 0 else 2, "solvent": "gas",
        "threads": 2, "psi4_memory_bytes": 512 * 1024**2,
        "timeout_seconds": 300, "process_tree_rss_limit_bytes": 1024**3,
        "sampling_interval_seconds": 0.1,
        "geometry_optimization": False, "frequency_calculation": False,
    }
    settings = root / "settings.json"
    settings.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    cmd = [args.psi4_python, str(root / "executed_psi4_single_point.py"), "--settings", str(settings)]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
    record = {
        "started_utc": datetime.now(timezone.utc).isoformat(), "command": cmd,
        "initial_available_memory_bytes": available, "config": cfg,
        "runner_sha256": sha256(__file__), "geometry_sha256": sha256(args.geometry),
        "status": "running", "memory_limit_kind": "sampled process-tree RSS, not a hard OS allocation quota",
    }
    started = time.perf_counter()
    peak = 0
    termination = None
    samples = []
    with (root / "stdout.log").open("w", encoding="utf-8") as stdout, (root / "stderr.log").open("w", encoding="utf-8") as stderr:
        proc = subprocess.Popen(cmd, cwd=root, stdout=stdout, stderr=stderr, env=env)
        parent = psutil.Process(proc.pid)
        while proc.poll() is None:
            elapsed = time.perf_counter() - started
            try:
                processes = [parent] + parent.children(recursive=True)
                rss = sum(p.memory_info().rss for p in processes if p.is_running())
            except psutil.NoSuchProcess:
                rss = 0
            peak = max(peak, rss)
            if not samples or elapsed - samples[-1][0] >= 1:
                samples.append([round(elapsed, 3), rss])
            if elapsed > cfg["timeout_seconds"]:
                termination = "timeout"
            elif rss > cfg["process_tree_rss_limit_bytes"]:
                termination = "memory_limit"
            if termination:
                try:
                    for child in parent.children(recursive=True):
                        child.kill()
                    parent.kill()
                except psutil.NoSuchProcess:
                    pass
                proc.wait(timeout=10)
                break
            time.sleep(cfg["sampling_interval_seconds"])
        exit_code = proc.wait()
    record.update(
        finished_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.perf_counter() - started, exit_code=exit_code,
        sampled_peak_process_tree_rss_bytes=peak, termination_reason=termination,
        status=termination or ("completed" if exit_code == 0 else "failed"),
        memory_samples_seconds_bytes=samples,
    )
    record["files"] = {
        str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*"))
        if p.is_file() and "scratch" not in p.parts
    }
    (root / "run_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: record[k] for k in ["status", "elapsed_seconds", "sampled_peak_process_tree_rss_bytes", "exit_code"]}))


if __name__ == "__main__":
    main()
