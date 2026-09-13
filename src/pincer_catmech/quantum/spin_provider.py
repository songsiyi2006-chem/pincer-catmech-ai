"""Matched-geometry spin-gradient providers with immutable native evidence.

Psi4 and PySCF are optional engines, imported only when selected. The Psi4
subprocess adapter permits a campaign runtime to use an existing separate
environment, with an enforceable timeout and no environment installation.
All public positions and returned analytic gradients use atomic units.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import uuid

import numpy as np

from .spin_mecp import SurfacePair


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def validate_spin_system(symbols, positions_bohr, charge, multiplicities, *, require_pair=True):
    # Avoid importing ASE into the separate Psi4 runtime.
    elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
    numbers = {symbol: i + 1 for i, symbol in enumerate(elements)}
    if not symbols or any(s not in numbers for s in symbols):
        raise ValueError("Unsupported or empty element list")
    x = np.asarray(positions_bohr, dtype=float)
    if x.shape != (len(symbols), 3) or not np.all(np.isfinite(x)):
        raise ValueError("Finite positions_bohr must match the element list")
    if isinstance(charge, bool) or not isinstance(charge, (int, np.integer)):
        raise ValueError("Integer charge required")
    electrons = sum(numbers[s] for s in symbols) - int(charge)
    if (electrons < 1 or not len(multiplicities) or len(set(multiplicities)) != len(multiplicities)
            or (require_pair and len(multiplicities) != 2)):
        raise ValueError("Distinct spin multiplicities and a positive electron count are required; pairs require two states")
    for multiplicity in multiplicities:
        if (isinstance(multiplicity, bool) or not isinstance(multiplicity, (int, np.integer))
                or multiplicity < 1 or multiplicity - 1 > electrons
                or (electrons - (multiplicity - 1)) % 2):
            raise ValueError("Multiplicity and electron count parity are inconsistent")
    if len(x) > 1:
        separation = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
        if np.any(separation[np.triu_indices(len(x), 1)] < 0.1):
            raise ValueError("Coincident nuclei are not permitted")
    return x


class Psi4SpinProvider:
    """Serial RKS/UKS (or RHF/UHF) analytic gradients at matched geometry.

    Psi4 has global state, so this in-process class must not be called from
    concurrent threads. Distinct scratch directories and fresh guesses avoid
    accidentally reusing the other spin state's orbitals.
    """
    def __init__(self, symbols, *, charge=0, multiplicities=(1, 3), method="pbe",
                 basis="def2-svp", workdir, memory_mib=500, threads=2,
                 scf_max_iterations=120, grid_radial=75, grid_spherical=194,
                 scf_algorithm="diis", soscf_start_convergence=1e-4,
                 e_convergence=1e-9, d_convergence=1e-8):
        if not 250 <= memory_mib <= 600 or not 1 <= threads <= 2:
            raise ValueError("Pilot provider requires 250-600 MiB configured memory and 1-2 threads")
        self.symbols, self.charge = list(symbols), charge
        self.multiplicities, self.method, self.basis = tuple(multiplicities), method, basis
        self.workdir = Path(workdir).resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.memory_mib, self.threads, self.scf_max_iterations = memory_mib, threads, scf_max_iterations
        if grid_radial < 20 or grid_spherical < 50 or scf_algorithm not in {"diis", "soscf"}:
            raise ValueError("Unsupported SCF algorithm or undersized diagnostic grid")
        self.grid_radial, self.grid_spherical, self.scf_algorithm = grid_radial, grid_spherical, scf_algorithm
        if not np.isfinite(soscf_start_convergence) or soscf_start_convergence <= 0:
            raise ValueError("Positive finite SOSCF switching residual required")
        self.soscf_start_convergence=soscf_start_convergence
        if not all(np.isfinite(v) and v>0 for v in (e_convergence,d_convergence)):
            raise ValueError("Positive finite SCF energy and density convergence tolerances required")
        self.e_convergence,self.d_convergence=e_convergence,d_convergence
        self.last_record = None

    def __call__(self, positions_bohr):
        x = validate_spin_system(self.symbols, positions_bohr, self.charge, self.multiplicities)
        record = self.evaluate_states(x)
        a, b = record["states"]
        return SurfacePair(a["energy_hartree"], b["energy_hartree"],
                           np.array(a["gradient_hartree_bohr"]), np.array(b["gradient_hartree_bohr"]),
                           f"multiplicity_{self.multiplicities[0]}", f"multiplicity_{self.multiplicities[1]}",
                           True, record)

    def evaluate_states(self, positions_bohr):
        """Evaluate one or more states serially; a state failure stops this call.

        The command-line single-state mode isolates failures when screening a
        matrix. Pair evaluation remains all-or-nothing for optimization.
        """
        x = validate_spin_system(self.symbols, positions_bohr, self.charge,
                                 self.multiplicities, require_pair=False)
        import psi4
        folder = self.workdir / f"pair_{uuid.uuid4().hex}"
        folder.mkdir(exist_ok=False)
        started = time.monotonic()
        request = dict(symbols=self.symbols, positions_bohr=x.tolist(), charge=self.charge,
                       multiplicities=self.multiplicities, method=self.method, basis=self.basis,
                       configured_memory_mib=self.memory_mib, threads=self.threads,
                       psi4_version=psi4.__version__, coordinate_units="bohr", energy_units="hartree",
                       gradient_units="hartree/bohr", scf_e_convergence=self.e_convergence,
                       scf_d_convergence=self.d_convergence,
                       scf_max_iterations=self.scf_max_iterations, grid_radial=self.grid_radial,
                       grid_spherical=self.grid_spherical, scf_algorithm=self.scf_algorithm,
                       soscf_start_convergence=self.soscf_start_convergence if self.scf_algorithm=="soscf" else None)
        geometry_bytes = np.asarray(x, dtype="<f8").tobytes() + " ".join(self.symbols).encode()
        request["geometry_sha256"] = hashlib.sha256(geometry_bytes).hexdigest()
        _write_json(folder / "request.json", request)
        psi4.set_num_threads(self.threads)
        psi4.set_memory(f"{self.memory_mib} MiB")
        state_records = []
        try:
            for multiplicity in self.multiplicities:
                state_dir = folder / f"multiplicity_{multiplicity}"
                state_dir.mkdir()
                psi4.core.clean()
                psi4.core.clean_options()
                psi4.core.clean_variables()
                psi4.core.IOManager.shared_object().set_default_path(str(state_dir) + os.sep)
                native_output = state_dir / "psi4.out"
                psi4.core.set_output_file(str(native_output), False)
                geometry = f"{self.charge} {multiplicity}\n" + "\n".join(
                    f"{symbol} {a:.16f} {b:.16f} {c:.16f}" for symbol, (a, b, c) in zip(self.symbols, x))
                geometry += "\nunits bohr\nno_reorient\nno_com\nsymmetry c1\n"
                mol = psi4.geometry(geometry)
                hf = self.method.lower() in {"hf", "scf"}
                reference = ("rhf" if hf else "rks") if multiplicity == 1 else ("uhf" if hf else "uks")
                psi4.set_options(dict(basis=self.basis, reference=reference, scf_type="df",
                                      guess="sad", e_convergence=self.e_convergence, d_convergence=self.d_convergence,
                                      maxiter=self.scf_max_iterations, fail_on_maxiter=True,
                                      dft_spherical_points=self.grid_spherical, dft_radial_points=self.grid_radial))
                if self.scf_algorithm == "soscf":
                    psi4.set_options(dict(soscf=True, soscf_start_convergence=self.soscf_start_convergence,
                                          soscf_max_iter=5, soscf_min_iter=1))
                state_start = time.monotonic()
                gradient, wfn = psi4.gradient(self.method, molecule=mol, return_wfn=True, dertype=1)
                g = np.asarray(gradient.to_array(), dtype=float)
                energy = float(wfn.energy())
                ca = np.asarray(wfn.Ca_subset("AO", "OCC").to_array())
                cb = np.asarray(wfn.Cb_subset("AO", "OCC").to_array())
                overlap = np.asarray(wfn.S().to_array())
                sz = (wfn.nalpha() - wfn.nbeta()) / 2
                s2 = float(sz * (sz + 1) + wfn.nbeta() - np.sum((ca.T @ overlap @ cb)**2))
                actual_x = np.asarray(mol.geometry().to_array())
                if not np.allclose(actual_x, x, atol=1e-10, rtol=0):
                    raise RuntimeError("Psi4 changed the geometry frame; matched gradients rejected")
                if not np.isfinite(energy) or not np.all(np.isfinite(g)):
                    raise RuntimeError("Psi4 returned nonfinite energy or gradient")
                psi4.core.flush_outfile()
                rec = dict(multiplicity=multiplicity, reference=reference, energy_hartree=energy,
                           gradient_hartree_bohr=g.tolist(), spin_squared=s2,
                           expected_spin_squared=(multiplicity-1)*(multiplicity+1)/4,
                           nalpha=wfn.nalpha(), nbeta=wfn.nbeta(), basis_functions=wfn.nmo(),
                           scf_converged=True, gradient_kind="analytic", wall_seconds=time.monotonic()-state_start,
                           output_sha256=hashlib.sha256(native_output.read_bytes()).hexdigest())
                _write_json(state_dir / "result.json", rec)
                state_records.append(rec)
            record = dict(**request, states=state_records, wall_seconds=time.monotonic()-started,
                          evidence="spin-sensitive electronic-structure calculation; no experimental validation",
                          native_directory=str(folder), status="completed")
            _write_json(folder / "result.json", record)
            self.last_record = record
            return record
        except Exception as exc:
            _write_json(folder / "failure.json", dict(error_type=type(exc).__name__, error=str(exc),
                         states_completed=state_records, wall_seconds=time.monotonic()-started))
            raise RuntimeError(f"Psi4 spin pair failed; inspect {folder}: {exc}") from exc
        finally:
            psi4.core.clean()


class Psi4SubprocessProvider:
    """Use a separate existing Psi4 Python, with per-pair deadline enforcement."""
    def __init__(self, python_executable, symbols, *, workdir, timeout_seconds=180, **settings):
        self.python_executable = str(Path(python_executable).resolve())
        if not Path(self.python_executable).is_file():
            raise FileNotFoundError(self.python_executable)
        if not np.isfinite(timeout_seconds) or timeout_seconds <= 0:
            raise ValueError("A positive finite timeout is required")
        self.symbols, self.settings = list(symbols), settings
        self.workdir, self.timeout_seconds = Path(workdir).resolve(), timeout_seconds
        self.workdir.mkdir(parents=True, exist_ok=True)

    def __call__(self, positions_bohr):
        x = validate_spin_system(self.symbols, positions_bohr, self.settings.get("charge", 0),
                                 self.settings.get("multiplicities", (1, 3)))
        folder = self.workdir / f"job_{uuid.uuid4().hex}"
        folder.mkdir(exist_ok=False)
        request = dict(symbols=self.symbols, positions_bohr=x.tolist(), **self.settings,
                       workdir=str(folder / "native"))
        _write_json(folder / "request.json", request)
        env = dict(os.environ)
        env.update(OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
                   PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
        # Conda DLL search must prefer the selected engine's environment.
        prefix = Path(self.python_executable).parent
        env["PATH"] = str(prefix / "Library" / "bin") + os.pathsep + str(prefix) + os.pathsep + env.get("PATH", "")
        env["PYTHONPATH"] = str(Path(__file__).resolve().parents[2])
        command = [self.python_executable, "-m", "pincer_catmech.quantum.spin_provider",
                   "--request", str(folder / "request.json"), "--result", str(folder / "response.json")]
        start = time.monotonic()
        with (folder / "worker.out").open("w", encoding="utf-8") as log:
            try:
                proc = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=self.timeout_seconds, cwd=folder)
            except subprocess.TimeoutExpired as exc:
                _write_json(folder / "launch.json", dict(argv=command, status="timeout", wall_seconds=time.monotonic()-start))
                raise TimeoutError(f"Psi4 pair exceeded {self.timeout_seconds} seconds; {folder}") from exc
        _write_json(folder / "launch.json", dict(argv=command, returncode=proc.returncode, wall_seconds=time.monotonic()-start))
        if proc.returncode or not (folder / "response.json").is_file():
            raise RuntimeError(f"Psi4 worker failed with exit {proc.returncode}; inspect {folder}")
        response = json.loads((folder / "response.json").read_text(encoding="utf-8"))
        return SurfacePair(**response).validate(x.shape)


class PySCFSpinProvider:
    """Optional analytic RKS/UKS or RHF/UHF implementation; no installation."""
    def __init__(self, symbols, *, charge=0, multiplicities=(1, 3), method="pbe",
                 basis="def2-svp", workdir, memory_mib=500, threads=2):
        self.symbols, self.charge, self.multiplicities = list(symbols), charge, tuple(multiplicities)
        self.method, self.basis, self.workdir = method, basis, Path(workdir)
        if not 250 <= memory_mib <= 600 or not 1 <= threads <= 2:
            raise ValueError("Pilot provider requires 250-600 MiB memory and 1-2 threads")
        self.memory_mib, self.threads = memory_mib, threads

    def __call__(self, positions_bohr):
        x = validate_spin_system(self.symbols, positions_bohr, self.charge, self.multiplicities)
        from pyscf import dft, gto, lib, scf
        lib.num_threads(self.threads)
        folder = self.workdir / f"pair_{uuid.uuid4().hex}"
        folder.mkdir(parents=True, exist_ok=False)
        records = []
        for multiplicity in self.multiplicities:
            mol = gto.M(atom=list(zip(self.symbols, x.tolist())), unit="Bohr", charge=self.charge,
                        spin=multiplicity-1, basis=self.basis, symmetry=False,
                        max_memory=self.memory_mib, verbose=4,
                        output=str(folder / f"spin_{multiplicity}.out"))
            if self.method.lower() in {"hf", "scf"}:
                mf = scf.RHF(mol) if multiplicity == 1 else scf.UHF(mol)
            else:
                mf = dft.RKS(mol) if multiplicity == 1 else dft.UKS(mol)
                mf.xc = self.method
                mf.grids.level = 3
            mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-10, 1e-7, 120
            energy = float(mf.kernel())
            if not mf.converged:
                raise RuntimeError(f"PySCF SCF not converged; {folder}")
            grad = np.asarray(mf.nuc_grad_method().kernel())
            s2 = 0.0 if multiplicity == 1 else float(mf.spin_square()[0])
            records.append(dict(multiplicity=multiplicity, energy_hartree=energy,
                                gradient_hartree_bohr=grad.tolist(), spin_squared=s2))
        metadata = dict(method=self.method, basis=self.basis, states=records, native_directory=str(folder),
                        symbols=self.symbols, positions_bohr=x.tolist(), charge=self.charge)
        _write_json(folder / "result.json", metadata)
        a, b = records
        return SurfacePair(a["energy_hartree"], b["energy_hartree"], np.array(a["gradient_hartree_bohr"]),
                           np.array(b["gradient_hartree_bohr"]), f"multiplicity_{self.multiplicities[0]}",
                           f"multiplicity_{self.multiplicities[1]}", True, metadata).validate(x.shape)


def _main():
    parser = argparse.ArgumentParser(description="Immutable single-process Psi4 spin-pair worker")
    parser.add_argument("--request", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()
    request = json.loads(Path(args.request).read_text(encoding="utf-8"))
    x = np.array(request.pop("positions_bohr"))
    single_state = request.pop("single_state", False)
    if single_state:
        if len(request["multiplicities"]) != 1:
            raise ValueError("single_state request requires exactly one multiplicity")
        _write_json(args.result, Psi4SpinProvider(**request).evaluate_states(x))
        return
    pair = Psi4SpinProvider(**request)(x)
    response = dict(energy1_hartree=pair.energy1_hartree, energy2_hartree=pair.energy2_hartree,
                    gradient1_hartree_bohr=pair.gradient1_hartree_bohr.tolist(),
                    gradient2_hartree_bohr=pair.gradient2_hartree_bohr.tolist(),
                    state1=pair.state1, state2=pair.state2, spin_dependent=pair.spin_dependent,
                    metadata=pair.metadata)
    _write_json(args.result, response)


if __name__ == "__main__":
    _main()
