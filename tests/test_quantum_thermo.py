"""Thermochemical identities and synthetic parser fixtures; no research numbers."""

import math
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms, units
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write

from pincer_catmech.kinetics.free_energy import concentration_correction, harmonic_vibrational_entropy
from pincer_catmech.quantum.thermochemistry import (
    DEFAULT_TEMPERATURES, EV_TO_KCAL, project_hessian, read_native_hessian,
    thermochemistry_from_modes,
)


def test_native_matrix_units_and_strict_shape(tmp_path):
    matrix = np.arange(36, dtype=float).reshape(6, 6) * .01
    path = tmp_path / "hessian.out"
    np.savetxt(path, matrix)
    readback = read_native_hessian(path, 2)
    assert readback == pytest.approx(matrix * units.Hartree / units.Bohr**2)
    path.write_text("1 2 3 NaN\n")
    with pytest.raises(ValueError):
        read_native_hessian(path, 2)


def test_project_unprojected_hessian_known_linear_bond():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    hessian = np.zeros((6, 6))
    hessian[2, 2] = hessian[5, 5] = 2
    hessian[2, 5] = hessian[5, 2] = -2
    result = project_hessian(atoms, hessian)
    assert result.external_rank == 5
    assert result.frequencies_cm1.shape == (1,)
    expected = math.sqrt(4 / atoms.get_masses()[0] * units._e / (units._amu * 1e-20)) / (2 * math.pi * units._c * 100)
    assert result.frequencies_cm1[0] == pytest.approx(expected)
    assert project_hessian(atoms, -hessian).frequencies_cm1[0] == pytest.approx(-expected)


def test_nearly_linear_geometry_uses_identical_projection_and_rrho_rank():
    atoms = Atoms("CO2", positions=[[0, 0, 0], [-1.2, 0, 0], [1.2, 1e-5, 0]])
    spectrum = project_hessian(atoms, np.eye(9))
    assert spectrum.external_rank == 6
    result = thermochemistry_from_modes(atoms, -10, spectrum.frequencies_cm1, temperatures=[298.15])[0]
    assert len(result.frequencies_cm1) == 3


def test_full_rrho_enthalpy_identity_and_explicit_temperature_grid():
    atoms = molecule("H2O")
    frequencies = np.array([1600., 3600., 3800.])
    energy = -10.0
    result = thermochemistry_from_modes(atoms, energy, frequencies, temperatures=[383.15])[0]
    vibration = frequencies * units.invcm
    thermal_vibration = np.sum(vibration / np.expm1(vibration / (units.kB * 383.15)))
    expected_h = energy + .5 * np.sum(vibration) + thermal_vibration + 4 * units.kB * 383.15
    assert result.H_298 == pytest.approx(expected_h * EV_TO_KCAL, rel=1e-12)
    assert result.temperature == 383.15
    assert result.S_harmonic > harmonic_vibrational_entropy(frequencies, 383.15)
    assert 383.15 in DEFAULT_TEMPERATURES and len(DEFAULT_TEMPERATURES) == 13
    assert not result.solvation_correction_added


def test_translational_rotational_and_electronic_entropy_against_independent_equations():
    atoms = molecule("H2O")
    temperature, pressure = 298.15, 101325.0
    frequencies = [1600, 3600, 3800]
    result = thermochemistry_from_modes(atoms, -10, frequencies, unpaired=2, symmetry_number=2, temperatures=[temperature])[0]
    mass = atoms.get_masses().sum() * units._amu
    moments = atoms.get_moments_of_inertia() * units._amu * 1e-20
    qtrans = (2 * math.pi * mass * units._k * temperature / units._hplanck**2)**1.5 * units._k * temperature / pressure
    qrot = math.sqrt(math.pi * np.prod(moments)) / 2 * (8 * math.pi**2 * units._k * temperature / units._hplanck**2)**1.5
    nonvib = units.kB * (math.log(qtrans) + 2.5 + math.log(qrot) + 1.5 + math.log(3)) * EV_TO_KCAL * 1000
    assert result.S_harmonic - result.S_vib_harmonic == pytest.approx(nonvib, abs=1e-5)


def test_one_molar_free_energy_pressure_invariance_and_no_double_shift():
    atoms = molecule("H2O")
    kwargs = dict(atoms=atoms, energy_ev=-10, frequencies_cm1=[1600, 3600, 3800], temperatures=[383.15])
    one = thermochemistry_from_modes(**kwargs, pressure_pa=101325)[0]
    two = thermochemistry_from_modes(**kwargs, pressure_pa=202650)[0]
    assert one.G_298_qRRHO_sol == pytest.approx(two.G_298_qRRHO_sol, abs=1e-6)
    assert one.G_298_qRRHO_sol - one.G_qRRHO_gas == pytest.approx(concentration_correction(383.15))
    with pytest.raises(ValueError, match="already include"):
        thermochemistry_from_modes(**kwargs, solvation_state="reference")


def test_ts_imaginary_mode_removed_once_without_silent_mode_selection():
    atoms = molecule("H2O")
    result = thermochemistry_from_modes(atoms, -10, [-600, 3600, 3800], transition_state=True, temperatures=[298.15])[0]
    assert result.frequencies_cm1 == (3600., 3800.)
    assert result.ZPVE == pytest.approx(.5 * (3600 + 3800) * units.invcm * EV_TO_KCAL)
    with pytest.raises(ValueError, match="all positive"):
        thermochemistry_from_modes(atoms, -10, [-5, 3600, 3800])
    with pytest.raises(ValueError, match="requires 3"):
        thermochemistry_from_modes(atoms, -10, [3600, 3800])
    with pytest.raises(ValueError, match="exactly one"):
        thermochemistry_from_modes(atoms, -10, [-600, -20, 3800], transition_state=True)


def test_explicit_non_dehydrogenation_ts_window_keeps_exact_negative_count():
    atoms = molecule("H2O")
    kwargs = dict(atoms=atoms, energy_ev=-10, frequencies_cm1=[-250, 3600, 3800], transition_state=True, temperatures=[383.15])
    with pytest.raises(ValueError, match="must lie"):
        thermochemistry_from_modes(**kwargs)
    result = thermochemistry_from_modes(**kwargs, ts_imaginary_window_cm1=None)[0]
    assert result.frequencies_cm1 == (3600., 3800.)
    assert any("One imaginary" in warning for warning in result.warnings)
    with pytest.raises(ValueError, match="exactly one"):
        thermochemistry_from_modes(atoms, -10, [-250, -5, 3800], transition_state=True, ts_imaginary_window_cm1=None)
    with pytest.raises(ValueError, match="ordered"):
        thermochemistry_from_modes(**kwargs, ts_imaginary_window_cm1=(-100, -1800))


@pytest.mark.parametrize("temperatures", [[], [298.15, 298.15], [-1], [math.inf]])
def test_invalid_temperature_grids(temperatures):
    with pytest.raises(ValueError):
        thermochemistry_from_modes(molecule("H2O"), -10, [1600, 3600, 3800], temperatures=temperatures)


def _stub_stationary_backend(monkeypatch, frequencies):
    """Explicit software doubles: no xTB call and no research dataset output."""
    import pincer_catmech.quantum.thermochemistry as bridge
    class ZeroForce(Calculator):
        implemented_properties = ["energy", "forces"]
        def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {"energy": -10., "forces": np.zeros((len(atoms), 3))}
    def mock_hessian(atoms, scratch_dir, **kwargs):
        folder = Path(scratch_dir)
        folder.mkdir()
        (folder / "hessian.out").write_text("synthetic archive fixture; never a research Hessian")
        (folder / "xtb.out").write_text("synthetic software fixture")
        (folder / "command.json").write_text(json.dumps({"test_double": True}))
        spectrum = project_hessian(atoms, np.eye(3 * len(atoms)))
        spectrum.frequencies_cm1 = np.asarray(frequencies, dtype=float)
        return spectrum, -10., 0.
    monkeypatch.setattr(bridge, "XTBCalculator", lambda *args, **kwargs: ZeroForce())
    monkeypatch.setattr(bridge, "native_hessian", mock_hessian)
    return bridge


def test_verified_minimum_source_hash_matches_referenced_spectrum(tmp_path, monkeypatch):
    bridge = _stub_stationary_backend(monkeypatch, [1600, 3600, 3800])
    result = bridge.characterize_stationary_point(molecule("H2O"), tmp_path / "output", tmp_path / "scratch",
                                                max_repairs=0, temperatures=[383.15])
    assert result["accepted"] and result["status"] == "verified_minimum"
    record = result["thermochemistry"][0]
    assert record["source_sha256"] == hashlib.sha256(Path(record["source"]).read_bytes()).hexdigest()
    assert record["source_sha256"] != result["attempts"][0]["unprojected_hessian"]["source_sha256"]


def test_changed_catalyst_identity_stops_before_hessian_or_thermochemistry(tmp_path, monkeypatch):
    import pincer_catmech.quantum.identity as identity
    bridge = _stub_stationary_backend(monkeypatch, [1600, 3600, 3800])
    monkeypatch.setattr(identity, "catalyst_identity_screen", lambda atoms, metadata: {"retained": False, "lost_bonds": [[0, 1]]})
    def forbidden(*args, **kwargs):
        raise AssertionError("Changed catalyst identity must not receive thermochemistry")
    monkeypatch.setattr(bridge, "native_hessian", forbidden)
    result = bridge.characterize_stationary_point(molecule("H2O"), tmp_path / "output", tmp_path / "scratch",
                                                identity_metadata={"test_double": True})
    assert not result["accepted"] and result["status"] == "chemical_identity_not_retained"
    assert result["thermochemistry"] == [] and len(result["attempts"]) == 1


def test_small_imaginary_modes_remain_rejected_after_exact_bounded_repairs(tmp_path, monkeypatch):
    bridge = _stub_stationary_backend(monkeypatch, [-5, 3600, 3800])
    calls = []
    @dataclass
    class Optimized:
        converged: bool = True
    def mock_optimize(atoms, folder, **kwargs):
        folder.mkdir()
        write(folder / "xtbopt.xyz", atoms)
        calls.append(kwargs)
        return Optimized()
    monkeypatch.setattr(bridge, "optimize_geometry", mock_optimize)
    result = bridge.characterize_stationary_point(molecule("H2O"), tmp_path / "output", tmp_path / "scratch",
                                                max_repairs=2, temperatures=[383.15])
    assert len(calls) == 2 and len(result["attempts"]) == 3
    assert not result["accepted"] and result["status"] == "numerical_uncertainty"
    assert result["thermochemistry"] == []
    assert all(attempt["imaginary_count"] == 1 and attempt["frequencies_cm1"][0] == -5 for attempt in result["attempts"])
    assert all(call["solvent"] == "toluene" and call["solvation_state"] == "gsolv" for call in calls)


@pytest.mark.parametrize("frozen,scale,valid", [(0, 1.0, True), (1, 1.0, False), (0, .95, False)])
def test_native_hessian_requires_full_unscaled_output_confirmation(tmp_path, monkeypatch, frozen, scale, valid):
    import pincer_catmech.quantum.thermochemistry as bridge
    def mock_execute(command, cwd, **kwargs):
        np.savetxt(Path(cwd) / "hessian.out", np.eye(6) * .01)
        return SimpleNamespace(returncode=0, stderr="", stdout=(
            f"Numerical Hessian\nfrozen atoms in % : {100*frozen/2:.5f} {frozen}\n"
            f"Hessian scale factor : {scale:.5f}\nTOTAL ENERGY -1.00000000 Eh\nnormal termination of xtb\n"))
    monkeypatch.setattr(bridge.subprocess, "run", mock_execute)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    def execute():
        return bridge.native_hessian(atoms, tmp_path / "native", executable=str(tmp_path / "test-double.exe"))
    if valid:
        spectrum, energy, _ = execute()
        assert spectrum.external_rank == 5 and energy == pytest.approx(-units.Hartree)
    else:
        with pytest.raises(ValueError, match="frozen|unscaled"):
            execute()
