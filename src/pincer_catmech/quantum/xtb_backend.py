"""ASE adapter and native optimization for the actual xTB executable.

GFN2-xTB energies are semiempirical, occupation-dependent and not reliable
ground-spin splittings. No surrogate or fabricated energy fallback is used.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write
from ase.units import Bohr, Hartree
import numpy as np

NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def executable_path(executable=None):
    selected = executable or os.environ.get("PINCER_XTB") or shutil.which("xtb")
    if not selected:
        raise FileNotFoundError("Provide PINCER_XTB or install xTB in the selected runtime")
    return str(Path(selected).resolve())


def calculation_environment(executable, threads):
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads),
               OPENBLAS_NUM_THREADS=str(threads), OMP_STACKSIZE="256M")
    env["PATH"] = str(Path(executable).parent) + os.pathsep + env.get("PATH", "")
    return env


def validate_electrons(atoms, charge, unpaired):
    if (isinstance(charge, bool) or isinstance(unpaired, bool) or
        not isinstance(charge, (int, np.integer)) or
        not isinstance(unpaired, (int, np.integer)) or unpaired < 0):
        raise ValueError("Charge and unpaired-electron count must be integers")
    electrons = int(sum(atoms.numbers)) - int(charge)
    if electrons <= 0 or unpaired > electrons or (electrons - unpaired) % 2:
        raise ValueError("Electron count and unpaired-electron parity are inconsistent")
    if np.any(atoms.pbc) or not np.all(np.isfinite(atoms.positions)):
        raise ValueError("An isolated finite geometry is required")


def parse_gradient(path, n_atoms):
    """Read the last Turbomole-style xTB gradient, returning eV/angstrom forces."""
    text = Path(path).read_text(encoding="utf-8")
    blocks = text.split("$grad")[1:]
    if not blocks:
        raise ValueError("xTB gradient file has no $grad block")
    body = blocks[-1].split("$end")[0].strip().splitlines()
    rows = body[-n_atoms:]
    if len(rows) != n_atoms:
        raise ValueError("Incomplete xTB gradient")
    gradient = np.array([[float(t.replace("D", "E").replace("d", "e")) for t in row.split()] for row in rows])
    if gradient.shape != (n_atoms, 3) or not np.all(np.isfinite(gradient)):
        raise ValueError("Malformed xTB gradient dimensions or values")
    return -gradient * Hartree / Bohr


def parse_energy(text):
    values = re.findall(rf"TOTAL ENERGY\s+({NUMBER})\s+Eh", text)
    if not values:
        raise ValueError("xTB did not print a final total energy")
    energy = float(values[-1].replace("D", "E").replace("d", "e")) * Hartree
    if not np.isfinite(energy):
        raise ValueError("Nonfinite xTB energy")
    return energy


def _execute(atoms, folder, executable, charge, unpaired, threads,
             arguments, timeout_seconds, electronic_temperature=300.0,
             solvent=None, solvation_state="gsolv", restart_file=None):
    validate_electrons(atoms, charge, unpaired)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    if any((folder / name).exists() for name in ("xtb.out", "gradient", "xtbopt.xyz")):
        raise FileExistsError("Use a fresh native calculation directory to preserve prior evidence")
    write(folder / "input.xyz", atoms, format="xyz")
    restart_metadata = None
    if restart_file is not None:
        restart_file = Path(restart_file)
        shutil.copy2(restart_file, folder / "xtbrestart")
        restart_metadata = {"source": str(restart_file), "sha256": hashlib.sha256(restart_file.read_bytes()).hexdigest()}
    # Only an explicitly provenance-tracked electronic restart is carried over.
    command = [executable, "input.xyz", "--gfn", "2", "--chrg", str(charge),
               "--uhf", str(unpaired), "--etemp", str(electronic_temperature),
               "--iterations", "500", "--acc", "0.5", *arguments]
    if solvent is not None:
        if solvation_state != "gsolv":
            raise ValueError("Use ALPB gsolv with the separate temperature-dependent 1atm->1M correction")
        command.extend(["--alpb", solvent, solvation_state])
    if restart_file is not None:
        command.append("--restart")
    atomic_json(folder / "launch.json", {"argv": command, "started_utc": utc_now(),
                "restart": restart_metadata, "threads": threads})
    start = time.monotonic()
    try:
        completed = subprocess.run(command, cwd=folder, env=calculation_environment(executable, threads),
                                   capture_output=True, text=True, encoding="utf-8", errors="replace",
                                   timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        def decode(value):
            return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""
        (folder / "xtb.out").write_text(decode(exc.stdout) + "\n" + decode(exc.stderr), encoding="utf-8")
        raise
    wall = time.monotonic() - start
    output = completed.stdout + "\n" + completed.stderr
    (folder / "xtb.out").write_text(output, encoding="utf-8")
    return completed.returncode, output, wall, command


class XTBCalculator(Calculator):
    """Real GFN2-xTB energies and analytic gradients exposed through ASE.

    Results always use 300 K electronic smearing. An elevated-smearing retry is
    diagnostic only and is followed by a new 300 K evaluation, never relabeled.
    Trace records retain every accepted geometry, energy and force array.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, workdir, *, charge=0, unpaired=0, threads=1,
                 executable=None, timeout_seconds=120, deadline=None, trace_path=None,
                 solvent=None, solvation_state="gsolv", **kwargs):
        super().__init__(**kwargs)
        self.workdir = Path(workdir)
        self.charge, self.unpaired, self.threads = charge, unpaired, threads
        self.executable = executable_path(executable)
        self.timeout_seconds, self.deadline = timeout_seconds, deadline
        self.solvent, self.solvation_state = solvent, solvation_state
        self.trace_path = Path(trace_path) if trace_path else self.workdir / "evaluations.jsonl.gz"
        self.evaluations = 0

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        if self.deadline is not None and time.time() >= self.deadline:
            raise TimeoutError("Campaign calculation deadline reached")
        attempts = []
        warm_restart = None
        for attempt, temperature in enumerate((300.0, 1000.0, 300.0)):
            folder = self.workdir / f"evaluation_{self.evaluations:06d}_{attempt}"
            allowed = self.timeout_seconds
            if self.deadline is not None:
                remaining = self.deadline - time.time()
                if remaining <= 0:
                    atomic_json(self.workdir / f"deadline_{self.evaluations:06d}.json", {"attempts": attempts})
                    raise TimeoutError("Campaign deadline exhausted during electronic retries")
                allowed = min(allowed, remaining)
            try:
                code, output, wall, command = _execute(self.atoms, folder, self.executable,
                    self.charge, self.unpaired, self.threads, ["--grad"], allowed, temperature,
                    self.solvent, self.solvation_state, restart_file=warm_restart)
                attempts.append({"attempt": attempt, "electronic_temperature_K": temperature,
                                 "returncode": code, "wall_seconds": wall,
                                 "output_sha256": hashlib.sha256(output.encode()).hexdigest()})
                if code or "normal termination of xtb" not in output.lower():
                    continue
                energy = parse_energy(output)
                forces = parse_gradient(folder / "gradient", len(self.atoms))
                if temperature != 300.0:
                    if (folder / "xtbrestart").is_file():
                        warm_restart = folder / "xtbrestart"
                    continue
                self.results = {"energy": energy, "forces": forces}
                self.trace_path.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(self.trace_path, "at", encoding="utf-8") as stream:
                    stream.write(json.dumps({"utc": utc_now(), "method": "GFN2-xTB",
                        "charge": self.charge, "unpaired_electrons": self.unpaired,
                        "solvent": self.solvent, "solvation_state": self.solvation_state,
                        "electronic_temperature_K": 300.0, "symbols": self.atoms.get_chemical_symbols(),
                        "positions_A": self.atoms.positions.tolist(), "energy_eV": energy,
                        "forces_eV_A": forces.tolist(), "attempts": attempts}, allow_nan=False) + "\n")
                return
            except (subprocess.TimeoutExpired, ValueError, OSError) as exc:
                attempts.append({"attempt": attempt, "error": str(exc)})
        atomic_json(self.workdir / f"failed_{self.evaluations:06d}.json", {"attempts": attempts})
        if self.deadline is not None and time.time() >= self.deadline:
            raise TimeoutError("Campaign deadline exhausted during real xTB evaluation")
        raise RuntimeError("Real xTB energy/gradient evaluation failed; inspect attempt records")


@dataclass
class OptimizationResult:
    status: str
    energy_eV: float | None
    max_force_eV_A: float | None
    wall_seconds: float
    output_directory: str
    charge: int
    unpaired_electrons: int
    converged: bool
    atoms_count: int
    output_sha256: str | None
    error: str | None = None


def optimize_geometry(atoms, output_dir, *, charge=0, unpaired=0, threads=1,
                      executable=None, max_cycles=250, timeout_seconds=600,
                      convergence="tight", restart_temperature=300.0,
                      solvent=None, solvation_state="gsolv", restart_file=None):
    """Run native GFN2-xTB optimization; retain all native files and failure status."""
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    if any((folder / name).exists() for name in ("xtb.out", "gradient", "xtbopt.xyz", "result.json")):
        raise FileExistsError("Use a fresh optimization directory; previous calculation records are immutable")
    executable = executable_path(executable)
    started = time.monotonic()
    try:
        code, output, wall, command = _execute(atoms, folder, executable, charge, unpaired,
            threads, ["--opt", convergence, "--cycles", str(max_cycles), "--grad"],
            timeout_seconds, restart_temperature, solvent, solvation_state, restart_file=restart_file)
        optimized = folder / "xtbopt.xyz"
        energy = parse_energy(output) if "TOTAL ENERGY" in output else None
        force = None
        if (folder / "gradient").is_file():
            force = float(np.linalg.norm(parse_gradient(folder / "gradient", len(atoms)), axis=1).max())
        converged = (code == 0 and optimized.is_file() and
                     "GEOMETRY OPTIMIZATION CONVERGED" in output and
                     "normal termination of xtb" in output.lower() and restart_temperature == 300.0)
        result = OptimizationResult("converged" if converged else "not_converged", energy, force,
            wall, str(folder), charge, unpaired, converged, len(atoms),
            hashlib.sha256(output.encode()).hexdigest())
        atomic_json(folder / "command.json", {"argv": command, "returncode": code,
                    "electronic_temperature_K": restart_temperature, "threads": threads,
                    "solvent": solvent, "solvation_state": solvation_state})
    except (subprocess.TimeoutExpired, ValueError, OSError) as exc:
        result = OptimizationResult("failed", None, None, time.monotonic()-started,
            str(folder), charge, unpaired, False, len(atoms), None, str(exc))
    atomic_json(folder / "result.json", asdict(result))
    return result
