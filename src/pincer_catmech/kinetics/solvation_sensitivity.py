"""Conditional cluster association sensitivity, never a bulk population model.

Uses existing full molecular RRHO enthalpies and spectra. It changes the qRRHO
entropy cutoff, not the electronic solvent model or the sampled conformers.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path, PurePosixPath

from .free_energy import R_J, calculate_thermochemistry

R_KCAL = R_J / 4184.0


def source_path(repo: Path, original: str) -> Path:
    """Relocate an archived source under data/, without trusting absolute paths."""
    parts = PurePosixPath(original.replace("\\", "/")).parts
    if "data" not in parts or ".." in parts:
        raise ValueError("An archive-relative data source is required")
    relative = Path(*parts[parts.index("data"):])
    result = (repo.resolve() / relative).resolve()
    if (repo.resolve() / "data") not in result.parents:
        raise ValueError("Source escaped the data archive")
    return result


def checked_source(repo: Path, record: dict) -> Path:
    path = source_path(repo, record["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError(f"Source hash mismatch: {path}")
    if "bytes" in record and path.stat().st_size != record["bytes"]:
        raise ValueError(f"Source size mismatch: {path}")
    return path


def recalculate_g(row: dict, cutoff_cm1: float) -> float:
    required_units = {"energy_unit": "kcal/mol", "entropy_unit": "cal/(mol K)",
                      "frequency_unit": "cm^-1", "temperature_unit": "K"}
    if any(row.get(k) != v for k, v in required_units.items()):
        raise ValueError("Thermochemistry units are absent or incompatible")
    if row["concentration_mol_l"] != 1.0:
        raise ValueError("This audit requires a dimensionless activity relative to 1 M")
    result = calculate_thermochemistry(
        row["frequencies_cm1"], E_elec=row["E_elec"], ZPVE=row["ZPVE"],
        H_298=row["H_298"], S_harmonic=row["S_harmonic"],
        temperature=row["temperature"], pressure_pa=row["pressure_pa"],
        concentration_mol_l=1.0, cutoff_cm1=cutoff_cm1,
        rotor_inertia=row["rotor_inertia"], source=row["source"],
        source_sha256=row["source_sha256"], program=row["program"],
        entropy_provenance=row["entropy_provenance"], warnings=tuple(row["warnings"]))
    return result.G_298_qRRHO_sol


def conditional_association(delta_g_standard: float, n_alcohol: int,
                            temperature: float, alcohol_activity: float) -> dict:
    """C + n A <=> CA_n; ln(a_CA_n/a_C) = -dG0/RT + n ln(a_A).

    This conditional activity ratio assumes equilibrium for this association
    only. It has no catalyst kinetics, no mass balance and no ensemble weights.
    """
    if (isinstance(n_alcohol, bool) or not isinstance(n_alcohol, int)
            or n_alcohol < 1):
        raise ValueError("Alcohol stoichiometry must be a positive integer")
    if not math.isfinite(delta_g_standard):
        raise ValueError("Finite free energy required")
    if any(not math.isfinite(x) or x <= 0 for x in (temperature, alcohol_activity)):
        raise ValueError("Temperature and dimensionless alcohol activity must be positive")
    effective = delta_g_standard - n_alcohol * R_KCAL * temperature * math.log(alcohol_activity)
    return {"effective_association_kcal_mol": effective,
            "log_activity_ratio_complex_to_bare": -effective / (R_KCAL * temperature)}
