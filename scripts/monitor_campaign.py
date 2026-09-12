"""Ten-minute local compute telemetry with same-trajectory energy drift."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json
import re
import time

import psutil

REPO = Path(__file__).resolve().parents[1]
TRACKED = ("run_conformer_campaign.py", "run_neb_campaign.py",
           "run_minimum_thermochemistry.py", "run_reference_thermochemistry.py",
           "run_condensation_search.py")


def snapshot():
    roots, processes = {}, {}
    for proc in psutil.process_iter(["pid", "name", "cmdline", "cpu_times", "memory_info"]):
        try:
            command = proc.info["cmdline"] or []
            script = next((x for x in TRACKED if any(x in arg for arg in command)), None)
            if not script or "python" not in (proc.info["name"] or "").lower():
                continue
            roots[proc.pid] = script
            for child in [proc, *proc.children(recursive=True)]:
                with child.oneshot():
                    processes[child.pid] = {"pid": child.pid, "task": script,
                        "name": child.name(), "cpu_seconds": sum(child.cpu_times()[:2]),
                        "rss_bytes": child.memory_info().rss}
        except (psutil.Error, OSError):
            continue
    data = REPO / "data" / "campaign_alpb_toluene"
    control = json.loads((data / "control.json").read_text())
    records = [json.loads(line) for line in (data / "records.jsonl").read_text().splitlines()]
    drift = []
    for record in reversed(records):
        if not record["native_result"]["converged"]:
            continue
        trajectory = Path(record["native_result"]["output_directory"]) / "xtbopt.log"
        if not trajectory.is_file():
            continue
        energies = re.findall(r"energy:\s*([-+0-9.Ee]+)", trajectory.read_text())
        if len(energies) > 1:
            drift.append({"catalyst_id": record["catalyst_id"], "state": record["state"],
                "conformer": record["conformer"], "last_iteration_delta_eV":
                (float(energies[-1]) - float(energies[-2])) * 27.211386245988,
                "native_max_force_eV_A": record["native_result"]["max_force_eV_A"],
                "definition": "Last minus penultimate energy from one native optimization trajectory"})
        if len(drift) >= 2:
            break
    return {"utc": datetime.now(timezone.utc).isoformat(), "controller": control,
            "same_optimization_energy_drift": drift,
            "local_compute_processes": list(processes.values()),
            "host_cpu_percent": psutil.cpu_percent(interval=0.25),
            "host_memory_available_bytes": psutil.virtual_memory().available,
            "hosted_subagent_CPU": "Unavailable: local process metrics above measure scientific workers only",
            "CPU_interpretation": "Cumulative lifetime per currently alive process; short-lived finished xTB jobs are not included"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deadline-epoch", required=True, type=float)
    parser.add_argument("--interval-seconds", default=600.0, type=float)
    args = parser.parse_args()
    if args.interval_seconds < 60:
        raise ValueError("telemetry interval must be at least 60 seconds")
    output = REPO / "data" / "campaign" / "combined_telemetry.jsonl"
    while True:
        start = time.time()
        try:
            record = snapshot()
        except (OSError, json.JSONDecodeError) as error:
            record = {"utc": datetime.now(timezone.utc).isoformat(), "status": "snapshot_retry_next_interval", "error": str(error)}
        with output.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
        print(json.dumps(record, allow_nan=False), flush=True)
        if time.time() >= args.deadline_epoch:
            break
        time.sleep(max(0.1, min(start + args.interval_seconds, args.deadline_epoch) - time.time()))


if __name__ == "__main__":
    main()
