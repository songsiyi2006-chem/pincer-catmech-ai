"""Numerical identities and defensive parsing on explicitly synthetic outputs."""

import math
from dataclasses import asdict

import pytest
from pydantic import TypeAdapter, ValidationError

from pincer_catmech.kinetics.free_energy import (
    B_AV, C_CM, H_PLANCK, HARTREE_TO_KCAL, K_B, R_CAL, R_J,
    ThermochemistryResult, calculate_thermochemistry, concentration_correction,
    harmonic_vibrational_entropy, parse_thermochemistry, qrrho_vibrational_entropy,
)

MODES = (25.4, 68.2, 120.5, 300.0, 500.0, 1000.0, 1600.0, 2900.0, 3100.0)


def gaussian_output(energy=-100.0, modes=MODES, temperature=298.15):
    """Format-realistic synthetic RRHO record, never a claimed quantum calculation."""
    lines = [" Entering Gaussian System, Link 0=g16", " Harmonic frequencies (cm**-1)"]
    for start in range(0, len(modes), 3):
        lines.append(" Frequencies -- " + " ".join(str(v) for v in modes[start:start + 3]))
    lines.extend([
        f" Temperature {temperature} Kelvin. Pressure 1.00000 Atm.",
        " Zero-point correction= 0.123456 (Hartree/Particle)",
        f" Sum of electronic and zero-point Energies= {energy + .123456:.6f}",
        f" Sum of electronic and thermal Enthalpies= {energy + .143456:.6f}",
        "                     E (Thermal)             CV                S",
        "                      KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin",
        " Total                  88.500             30.000            120.000",
        " Electronic              0.000              0.000              0.000",
        " Normal termination of Gaussian 16 at Sat Sep 12 00:00:00 2026.",
    ])
    return "\n".join(lines) + "\n"


def orca_output(modes=MODES, external=6, builtin_vib_ts=.01000000):
    lines = [" * O   R   C   A *", " VIBRATIONAL FREQUENCIES", " -----------------------"]
    for i, value in enumerate((0.0,) * external + tuple(modes)):
        lines.append(f" {i:3d}: {value:12.6f} cm**-1")
    lines.extend([
        " THERMOCHEMISTRY AT 298.15K", " Temperature ... 298.15 K",
        " Pressure ... 1.00 atm", " Electronic energy ... -100.00000000 Eh",
        " Zero point energy ... 0.12345600 Eh 77.47 kcal/mol",
        " Total Enthalpy ... -99.85654400 Eh", " ENTROPY",
        " Vibrational entropy computed according to the QRRHO of S. Grimme",
        " Electronic entropy ... 0.00010000 Eh 0.06 kcal/mol",
        f" Vibrational entropy ... {builtin_vib_ts:.8f} Eh 6.28 kcal/mol",
        " Rotational entropy ... 0.00800000 Eh 5.02 kcal/mol",
        " Translational entropy ... 0.01800000 Eh 11.30 kcal/mol",
        " Final entropy term ... 0.03610000 Eh 22.65 kcal/mol",
        " GIBBS FREE ENERGY", " Total enthalpy ... -99.85654400 Eh",
        " ****ORCA TERMINATED NORMALLY****",
    ])
    return "\n".join(lines) + "\n"


def manual(**kwargs):
    values = dict(frequencies_cm1=MODES, E_elec=-1000.0, ZPVE=15.0,
                  H_298=-982.0, S_harmonic=120.0)
    values.update(kwargs)
    return calculate_thermochemistry(**values)


def test_independent_mode_equations_and_reduced_inertia():
    frequency, temperature = 68.2, 298.15
    nu_hz = frequency * C_CM
    x = H_PLANCK * nu_hz / (K_B * temperature)
    s_ho = R_CAL * (x / (math.exp(x) - 1) - math.log(1 - math.exp(-x)))
    mu = H_PLANCK / (8 * math.pi**2 * nu_hz)
    mu_eff = mu * B_AV / (mu + B_AV)
    s_rot = R_CAL * (0.5 + math.log(math.sqrt(8 * math.pi**3 * mu_eff * K_B * temperature / H_PLANCK**2)))
    weight = 1 / (1 + (100 / frequency)**4)
    assert harmonic_vibrational_entropy([frequency]) == pytest.approx(s_ho, abs=1e-12)
    assert qrrho_vibrational_entropy([frequency]) == pytest.approx(weight * s_ho + (1 - weight) * s_rot, abs=1e-12)


def test_full_entropy_and_enthalpy_preserved():
    result = manual()
    assert result.S_qRRHO < result.S_harmonic
    assert result.S_harmonic - result.S_vib_harmonic == pytest.approx(result.S_qRRHO - result.S_vib_qRRHO)
    assert result.H_298 == -982.0
    assert result.G_298_qRRHO_sol == pytest.approx(result.H_298 - result.temperature * result.S_qRRHO / 1000 + result.delta_G_concentration)
    assert not result.solvation_correction_added
    assert TypeAdapter(ThermochemistryResult).validate_json(TypeAdapter(ThermochemistryResult).dump_json(result)) == result


def test_exact_concentration_shift_and_pressure_scaling():
    expected = R_J * 298.15 / 4184 * math.log(1000 * R_J * 298.15 / 101325)
    assert concentration_correction() == pytest.approx(expected, abs=1e-14)
    assert 1.894 < expected < 1.895
    assert concentration_correction(298.15, 202650) == pytest.approx(expected - R_J * 298.15 / 4184 * math.log(2))
    assert concentration_correction(310, 101325, .01) == pytest.approx(R_J * 310 / 4184 * math.log(10 * R_J * 310 / 101325))


def test_entropy_stability_extreme_modes_and_grimme_low_frequency_limit():
    assert math.isfinite(harmonic_vibrational_entropy([1e-200, 1e200]))
    assert harmonic_vibrational_entropy([1e200]) == 0.0
    capped = qrrho_vibrational_entropy([1e-200])
    limit = R_CAL * (.5 + .5 * math.log(8 * math.pi**3 * B_AV * K_B * 298.15 / H_PLANCK**2))
    assert capped == pytest.approx(limit, abs=1e-12)
    assert qrrho_vibrational_entropy([1e-200], rotor_inertia="uncapped") > capped
    assert "Uncapped" in " ".join(manual(rotor_inertia="uncapped").warnings)


@pytest.mark.parametrize("kwargs", [
    {"temperature": 0}, {"pressure_pa": -1}, {"concentration_mol_l": math.inf},
    {"cutoff_cm1": 0}, {"frequencies_cm1": []}, {"frequencies_cm1": [math.nan]},
    {"ZPVE": -1}, {"H_298": math.inf}, {"S_harmonic": 1}, {"rotor_inertia": "unknown"},
])
def test_invalid_inputs_fail_explicitly(kwargs):
    with pytest.raises((ValueError, ValidationError)):
        manual(**kwargs)


def test_imaginary_and_zero_modes_require_deliberate_policy():
    with pytest.raises(ValueError, match="imaginary"):
        manual(frequencies_cm1=(-250.0,) + MODES)
    transition_state = manual(frequencies_cm1=(-250.0,) + MODES, transition_state=True)
    assert transition_state.frequencies_cm1 == MODES
    assert transition_state.S_qRRHO == manual().S_qRRHO
    with pytest.raises(ValueError, match="exactly one"):
        manual(transition_state=True)
    with pytest.raises(ValueError, match="exactly one"):
        manual(frequencies_cm1=(-50, -250) + MODES, transition_state=True)
    with pytest.raises(ValueError, match="Zero frequencies"):
        manual(frequencies_cm1=(0,) + MODES)
    assert manual(frequencies_cm1=(0,) + MODES, exclude_zero_modes=True).frequencies_cm1 == MODES


def test_result_schema_rejects_nonfinite_data():
    data = asdict(manual())
    data["G_298_qRRHO_sol"] = math.nan
    with pytest.raises(ValidationError):
        ThermochemistryResult(**data)


def test_gaussian_direct_parse(tmp_path):
    path = tmp_path / "ru_pnp_mock.log"
    path.write_text(gaussian_output(), encoding="utf-8")
    result = parse_thermochemistry(path)
    assert result.E_elec == pytest.approx(-100 * HARTREE_TO_KCAL, abs=1e-8)
    assert result.ZPVE == pytest.approx(.123456 * HARTREE_TO_KCAL)
    assert result.S_harmonic == 120
    assert result.program == "Gaussian"
    assert len(result.source_sha256) == 64
    assert result.job_index == 0


def test_orca_reconstructs_harmonic_entropy_without_double_correction(tmp_path):
    path = tmp_path / "ru_pnp_mock.out"
    path.write_text(orca_output(), encoding="utf-8")
    result = parse_thermochemistry(path)
    nonvib = .0261 * HARTREE_TO_KCAL * 1000 / 298.15
    assert result.S_harmonic == pytest.approx(nonvib + harmonic_vibrational_entropy(MODES))
    assert result.S_qRRHO == pytest.approx(nonvib + qrrho_vibrational_entropy(MODES))
    assert "6 identified zero external" in " ".join(result.warnings)
    path.write_text(orca_output(builtin_vib_ts=.999), encoding="utf-8")
    assert parse_thermochemistry(path).S_qRRHO == result.S_qRRHO


def test_orca_missing_entropy_components_and_unidentified_zeros_rejected(tmp_path):
    path = tmp_path / "bad.out"
    path.write_text(orca_output().replace(" Translational entropy", " Missing entropy"), encoding="utf-8")
    with pytest.raises(ValueError, match="reconstruct harmonic entropy"):
        parse_thermochemistry(path)
    path.write_text(orca_output(external=1), encoding="utf-8")
    with pytest.raises(ValueError, match="projected external"):
        parse_thermochemistry(path)


@pytest.mark.parametrize("factory,suffix,termination", [
    (gaussian_output, ".log", " Normal termination of Gaussian"),
    (orca_output, ".out", " ****ORCA TERMINATED NORMALLY"),
])
def test_truncated_termination_and_incomplete_final_jobs(tmp_path, factory, suffix, termination):
    path = tmp_path / ("incomplete" + suffix)
    truncated = factory().split(termination)[0]
    path.write_text(truncated, encoding="utf-8")
    with pytest.raises(ValueError, match="normal termination"):
        parse_thermochemistry(path)
    path.write_text(factory() + truncated, encoding="utf-8")
    with pytest.raises(ValueError, match="final job is incomplete"):
        parse_thermochemistry(path)
    assert parse_thermochemistry(path, job_index=0).job_index == 0


def test_gaussian_multijob_selection_never_mixes_records(tmp_path):
    path = tmp_path / "multiple.log"
    path.write_text(gaussian_output(-100) + gaussian_output(-200), encoding="utf-8")
    assert parse_thermochemistry(path).E_elec == pytest.approx(-200 * HARTREE_TO_KCAL)
    assert parse_thermochemistry(path, job_index=0).E_elec == pytest.approx(-100 * HARTREE_TO_KCAL)
    broken_final = gaussian_output(-200).replace(" Zero-point correction=", " Missing correction=")
    path.write_text(gaussian_output(-100) + broken_final, encoding="utf-8")
    with pytest.raises(ValueError, match="ZPVE"):
        parse_thermochemistry(path)


def test_temperature_mismatch_and_fortran_exponents(tmp_path):
    path = tmp_path / "exponents.log"
    path.write_text(gaussian_output().replace("0.123456", "1.23456D-01"), encoding="utf-8")
    assert parse_thermochemistry(path, temperature=298.15).ZPVE == pytest.approx(.123456 * HARTREE_TO_KCAL)
    with pytest.raises(ValueError, match="Requested temperature"):
        parse_thermochemistry(path, temperature=310)


def test_orca_multiple_temperature_tables_are_not_silently_combined(tmp_path):
    path = tmp_path / "multi-temperature.out"
    path.write_text(orca_output().replace(" GIBBS FREE ENERGY", " THERMOCHEMISTRY AT 310.00K\n GIBBS FREE ENERGY"), encoding="utf-8")
    with pytest.raises(ValueError, match="one thermochemistry temperature"):
        parse_thermochemistry(path)


def test_source_mode_filtering_is_not_silently_mixed(tmp_path):
    path = tmp_path / "filtered.out"
    output = orca_output().replace(" VIBRATIONAL FREQUENCIES", " CutOffFreq 35\n VIBRATIONAL FREQUENCIES")
    path.write_text(output, encoding="utf-8")
    with pytest.raises(ValueError, match="CutOffFreq"):
        parse_thermochemistry(path)
    path.write_text(orca_output(modes=(.5,) + MODES), encoding="utf-8")
    with pytest.raises(ValueError, match="CutOffFreq"):
        parse_thermochemistry(path)
    output = orca_output().replace(" Electronic energy", " freq. 68.20 E(vib) ... 0.10\n Electronic energy")
    path.write_text(output, encoding="utf-8")
    with pytest.raises(ValueError, match="thermal vibration list"):
        parse_thermochemistry(path)


def test_gaussian_printed_vibrational_entropy_checks_scaling(tmp_path):
    path = tmp_path / "scaled.log"
    output = gaussian_output().replace(" Electronic", " Vibrational            10.0              10.0             99.0\n Electronic")
    path.write_text(output, encoding="utf-8")
    with pytest.raises(ValueError, match="vibrational entropy disagrees"):
        parse_thermochemistry(path)


def test_incomplete_preceding_job_cannot_supply_final_job_thermochemistry(tmp_path):
    path = tmp_path / "concatenated.log"
    failed = gaussian_output().split(" Normal termination")[0] + " Error termination\n"
    final_sp = " Entering Gaussian System, Link 0=g16\n SCF Done: E(RB3LYP) = -200.0\n Normal termination of Gaussian 16\n"
    path.write_text(failed + final_sp, encoding="utf-8")
    with pytest.raises(ValueError, match="frequency block"):
        parse_thermochemistry(path)


@pytest.mark.parametrize("damaged", ["******", "68.2D+", "NaN", "inf"])
def test_gaussian_damaged_mode_tokens_are_not_skipped(tmp_path, damaged):
    path = tmp_path / "damaged.log"
    path.write_text(gaussian_output().replace("68.2", damaged), encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed Gaussian frequency"):
        parse_thermochemistry(path)


def test_orca_damaged_final_mode_is_not_skipped(tmp_path):
    path = tmp_path / "damaged.out"
    path.write_text(orca_output().replace("3100.000000", "***********"), encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed ORCA indexed frequency"):
        parse_thermochemistry(path)


def test_gaussian_empty_frequency_row_is_not_skipped(tmp_path):
    path = tmp_path / "empty-row.log"
    path.write_text(gaussian_output().replace(" Frequencies -- 1600.0 2900.0 3100.0", " Frequencies -- "), encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed Gaussian frequency"):
        parse_thermochemistry(path)
