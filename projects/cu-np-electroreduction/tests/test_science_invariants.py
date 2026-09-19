"""SYNTHETIC SOFTWARE TESTS of scientific conservation, not research results."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from science import KB_EV, F_C_MOL, rate_pair, site_steady_state, electrical_kwh_kg, convert_potential


class IndependentInvariantChecks(unittest.TestCase):
    def test_nonequal_site_energies_give_boltzmann_equilibrium(self):
        states = [0.0, 0.03, 0.05]
        rates = (*rate_pair(states[0], states[1], 0.20),
                 *rate_pair(states[1], states[2], 0.31),
                 *rate_pair(states[2], states[0], 0.40))
        populations, flux = site_steady_state(rates)
        weights = [math.exp(-g / (KB_EV * 298.15)) for g in states]
        expected = [x / sum(weights) for x in weights]
        for actual, reference in zip(populations, expected):
            self.assertAlmostEqual(actual, reference, places=11)
        self.assertLess(abs(flux), 1e-12 * max(rates))

    def test_flux_is_identical_through_all_three_edges(self):
        rates = [0.5, 0.01, 3.0, 0.8, 0.06, 0.03]
        p, flux = site_steady_state(rates)
        self.assertAlmostEqual(flux, rates[2] * p[1] - rates[3] * p[2], places=12)
        self.assertAlmostEqual(flux, rates[4] * p[2] - rates[5] * p[0], places=12)

    def test_full_cell_energy_from_explicit_one_kilogram_electron_balance(self):
        # Exactly 1 kg product, M=250 g/mol, z=4, FE=1/2, recovery=4/5.
        acceptable_moles = 1000.0 / 250.0
        target_formed_moles = acceptable_moles / 0.8
        electron_moles = target_formed_moles * 4.0 / 0.5
        electrical_joules = electron_moles * F_C_MOL * 5.0
        self.assertAlmostEqual(electrical_kwh_kg(5.0, 4.0, 0.5, 250.0, 0.8), electrical_joules / 3.6e6, places=12)

    def test_potential_reference_conversion_is_round_trip_consistent(self):
        a = dict(solvent="DMF", temperature_K=298.15, common_reference="X", offset_V=0.12, source="synthetic measured calibration")
        b = {**a, "offset_V": -0.31}
        shifted = convert_potential(-1.2, a, b)
        self.assertAlmostEqual(convert_potential(shifted, b, a), -1.2, places=12)
        self.assertAlmostEqual(shifted + b["offset_V"], -1.2 + a["offset_V"], places=12)


if __name__ == "__main__":
    unittest.main()
