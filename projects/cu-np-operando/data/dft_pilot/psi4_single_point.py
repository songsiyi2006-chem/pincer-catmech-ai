"""One fixed-nuclei gas-phase DFT diagnostic; run by bounded_runner.py."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import traceback

import numpy as np
import psi4


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--settings", required=True, type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.settings.read_text(encoding="utf-8"))
    root = args.settings.resolve().parent
    xyz = root / "input.xyz"
    psi4.set_num_threads(cfg["threads"])
    psi4.set_memory(cfg["psi4_memory_bytes"])
    psi4.core.set_output_file(str(root / "psi4.out"), False)
    scratch = root / "scratch"
    scratch.mkdir(exist_ok=True)
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    geom = f'{cfg["charge"]} {cfg["multiplicity"]}\n' + "\n".join(xyz.read_text().splitlines()[2:])
    geom += "\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"
    mol = psi4.geometry(geom)
    options = {
        "basis": cfg["basis"], "scf_type": "df",
        "df_basis_scf": "def2-universal-jkfit",
        "reference": "rks" if cfg["multiplicity"] == 1 else "uks",
        "e_convergence": 1e-8, "d_convergence": 1e-8,
        "maxiter": 150, "guess": "sad", "fail_on_maxiter": True,
        "dft_radial_points": 75, "dft_spherical_points": 302,
        "ints_tolerance": 1e-12,
    }
    if cfg.get("soscf_rescue", False):
        options.update(soscf=True, soscf_start_convergence=1e-3, soscf_max_iter=5)
    psi4.set_options(options)
    result = {
        "evidence_class": "CALCULATION", "status": "started",
        "scope": "gas-phase fixed-nuclei molecular diagnostic, no optimization or frequency calculation",
        "settings": cfg, "psi4_version": psi4.__version__,
        "input_sha256": digest(xyz), "settings_sha256": digest(args.settings),
        "child_script_sha256": digest(__file__), "psi4_options": options,
    }
    started = time.perf_counter()
    try:
        energy, wfn = psi4.energy(cfg["method"], molecule=mol, return_wfn=True)
        na, nb = wfn.nalpha(), wfn.nbeta()
        ca = wfn.Ca().np[:, :na]
        cb = wfn.Cb().np[:, :nb]
        overlap = wfn.S().np
        spin_z = (na - nb) / 2.0
        s2 = spin_z * (spin_z + 1.0) + nb - np.sum((ca.T @ overlap @ cb) ** 2)
        alpha = wfn.epsilon_a().np.tolist()
        beta = wfn.epsilon_b().np.tolist()
        basis_path = Path(psi4.core.get_datadir()) / "basis" / (cfg["basis"].lower() + ".gbs")
        result.update({
            "status": "completed", "energy_hartree": float(energy),
            "nalpha": na, "nbeta": nb, "nbf": wfn.basisset().nbf(),
            "s2_orbital_overlap": float(s2),
            "s2_ideal": spin_z * (spin_z + 1.0),
            "s2_formula": "Sz*(Sz+1)+Nbeta-sum(abs(Ca_occ.T*S*Cb_occ)**2), unrestricted determinant diagnostic",
            "orbital_energies_alpha_hartree": alpha,
            "orbital_energies_beta_hartree": beta,
            "alpha_homo_hartree": alpha[na - 1], "alpha_lumo_hartree": alpha[na],
            "beta_homo_hartree": beta[nb - 1], "beta_lumo_hartree": beta[nb],
            "basis_source_sha256": digest(basis_path),
            "diffuse_basis": "def2-SVPD includes diffuse augmentation; no diffuse-basis convergence claim",
        })
        try:
            psi4.molden(wfn, str(root / "orbitals.molden"))
            result["molden_written"] = True
        except Exception as exc:
            result["molden_written"] = False
            result["molden_error"] = repr(exc)
    except Exception as exc:
        result.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        (root / "python_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise
    finally:
        result["child_elapsed_seconds"] = time.perf_counter() - started
        (root / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        psi4.core.clean()


if __name__ == "__main__":
    main()
