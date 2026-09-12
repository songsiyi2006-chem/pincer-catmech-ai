"""End-to-end synthetic Ru-PNP-like fixtures, never experimental evidence."""

from dataclasses import asdict
import json
import math

from ase import Atoms
import numpy as np
from pydantic import TypeAdapter
import pytest

from pincer_catmech.features.steric import bite_angle, buried_volume
from pincer_catmech.kinetics.free_energy import (
    ThermochemistryResult, calculate_thermochemistry, concentration_correction,
    harmonic_vibrational_entropy, parse_thermochemistry,
)


@pytest.fixture
def ru_pnp():
    """Toy scaffold and 3N-6 assigned modes; no claim of chemical connectivity."""
    atoms = Atoms(
        ["Ru", "P", "P", "N", "C", "C", "H", "H"], positions=[
            [0, 0, 0], [2.3, 0, 0], [-2.3, 0, 0], [0, 2.1, 0],
            [1.6, 1.5, 0], [-1.6, 1.5, 0], [1.7, 1.5, 1.1], [-1.7, 1.5, 1.1],
        ],
    )
    modes = [25.4, 68.2, 120.5, 160, 210, 300, 420, 550, 700,
             850, 1000, 1100, 1250, 1400, 1500, 1650, 3000, 3100]
    return atoms, modes


def test_mock_ru_pnp_complete_pipeline(ru_pnp):
    atoms, modes = ru_pnp
    assert len(modes) == 3 * len(atoms) - 6
    entropy = 65.0 + harmonic_vibrational_entropy(modes)
    thermo = calculate_thermochemistry(
        modes, E_elec=-1000000, ZPVE=25, H_298=-999970,
        S_harmonic=entropy, source="synthetic Ru-PNP fixture",
    )
    assert thermo.S_qRRHO < thermo.S_harmonic
    assert thermo.G_298_qRRHO_sol == pytest.approx(
        -999970 - 298.15 * thermo.S_qRRHO / 1000 + concentration_correction()
    )
    assert TypeAdapter(ThermochemistryResult).validate_json(
        TypeAdapter(ThermochemistryResult).dump_json(thermo)
    ) == thermo
    assert bite_angle(atoms, 0, 1, 2) == pytest.approx(180.0)
    assert bite_angle(atoms, 0, 1, 3) == pytest.approx(90.0)
    feature = buried_volume(atoms, 0, ligand_indices=range(1, len(atoms)), seed=2026)
    low, high = feature.confidence_interval_percent
    assert 0 < low < feature.percent_buried_volume < high < 100
    payload = {"thermo": asdict(thermo), "steric": asdict(feature)}
    assert json.loads(json.dumps(payload, allow_nan=False))["thermo"]["program"] == "manual"


def test_analytical_oblique_bite_angle():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, math.sqrt(3), 0.0]])
    assert bite_angle(coords, 0, 1, 2, symbols=["Ru", "P", "N"]) == pytest.approx(60.0)


def test_gaussian_mock_file_to_features(tmp_path, ru_pnp):
    atoms, modes = ru_pnp
    entropy = 65.0 + harmonic_vibrational_entropy(modes)
    frequency_lines = "\n".join(
        " Frequencies -- " + " ".join(str(value) for value in modes[i:i+3])
        for i in range(0, len(modes), 3)
    )
    source = tmp_path / "mock_ru_pnp.log"
    source.write_text(
        "Entering Gaussian System, Link 0=g16\n"
        "SYNTHETIC DATA: not an executed Gaussian job\n"
        f"{frequency_lines}\n"
        " Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.\n"
        " Zero-point correction=                           0.040000 (Hartree/Particle)\n"
        " Sum of electronic and zero-point Energies=       -1593.560000\n"
        " Sum of electronic and thermal Enthalpies=        -1593.550000\n"
        "                     E (Thermal)             CV                S\n"
        "                      KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin\n"
        f" Total                   32.0                40.0             {entropy:.9f}\n"
        " Normal termination of Gaussian 16\n", encoding="utf-8",
    )
    parsed = parse_thermochemistry(source)
    assert parsed.S_qRRHO < parsed.S_harmonic
    assert parsed.frequencies_cm1 == tuple(modes)
    assert len(parsed.source_sha256) == 64
    assert bite_angle(atoms, 0, 1, 3) == pytest.approx(90.0)
    assert math.isfinite(parsed.G_298_qRRHO_sol)
