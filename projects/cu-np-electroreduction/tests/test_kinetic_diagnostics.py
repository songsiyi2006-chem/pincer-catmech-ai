"""SYNTHETIC SOFTWARE TESTS ONLY; no Cu reaction labels or observations."""

import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from kinetic_diagnostics import (
    SYNTHETIC_LABEL, evaluate_network, perturb_transition_state,
    transition_state_rate_control, propagate_barrier_uncertainty,
)
from science import KB_EV, rate_pair


class KineticDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.rates = (0.5, 0.01, 3.0, 0.8, 0.06, 0.03)
        self.label = {"evidence_label": SYNTHETIC_LABEL}

    def test_common_transition_state_shift_preserves_ratios_and_affinity(self):
        for edge in range(3):
            shifted = perturb_transition_state(self.rates, edge, 0.04, **self.label)
            for i in range(3):
                self.assertAlmostEqual(shifted[2*i] / shifted[2*i+1], self.rates[2*i] / self.rates[2*i+1], places=12)
            before = evaluate_network(self.rates, **self.label)
            after = evaluate_network(shifted, **self.label)
            self.assertAlmostEqual(before["cycle_affinity_over_kT"], after["cycle_affinity_over_kT"], places=12)
            self.assertLess(shifted[2*edge], self.rates[2*edge])

    def test_equilibrium_flux_zero_and_boltzmann_populations(self):
        # Independent thermodynamic construction with unequal state energies.
        energies = (0.0, 0.03, 0.05)
        equilibrium = (*rate_pair(energies[0], energies[1], 0.20),
                       *rate_pair(energies[1], energies[2], 0.31),
                       *rate_pair(energies[2], energies[0], 0.40))
        result = propagate_barrier_uncertainty(
            equilibrium, (0.01, 0.03, 0.02), samples=30, seed=19,
            shared_half_width_eV=0.01, **self.label
        )
        weights = [math.exp(-energy / (KB_EV * 298.15)) for energy in energies]
        expected = [x / sum(weights) for x in weights]
        for draw in result["draws"]:
            self.assertLess(abs(draw["flux_per_site_s"]), 1e-12 * max(draw["rates_per_s"]))
            self.assertLess(draw["relative_net_flux"], 1e-12)
            for actual, reference in zip(draw["populations"], expected):
                self.assertAlmostEqual(actual, reference, places=11)
        with self.assertRaisesRegex(ValueError, "equilibrium"):
            transition_state_rate_control(equilibrium, **self.label)

    def test_site_balance_and_equal_edge_fluxes(self):
        result = propagate_barrier_uncertainty(
            self.rates, (0.03, 0.02, 0.04), samples=50, seed=11, **self.label
        )
        for draw in result["draws"]:
            self.assertAlmostEqual(sum(draw["populations"]), 1.0, places=13)
            self.assertGreaterEqual(min(draw["populations"]), 0.0)
            for flux in draw["edge_fluxes_per_site_s"]:
                self.assertAlmostEqual(flux, draw["flux_per_site_s"], places=12)
            for i, ratio in enumerate(draw["log_forward_reverse_ratios"]):
                self.assertAlmostEqual(ratio, math.log(self.rates[2*i] / self.rates[2*i+1]), places=12)

    def test_fixed_seed_reproduces_all_draws_and_different_seed_changes(self):
        kwargs = dict(samples=25, seed=731, shared_half_width_eV=0.01, **self.label)
        first = propagate_barrier_uncertainty(self.rates, (0.01, 0.02, 0.03), **kwargs)
        second = propagate_barrier_uncertainty(self.rates, (0.01, 0.02, 0.03), **kwargs)
        self.assertEqual(first, second)
        kwargs["seed"] = 732
        other = propagate_barrier_uncertainty(self.rates, (0.01, 0.02, 0.03), **kwargs)
        self.assertNotEqual(first["draws"], other["draws"])
        for draw in first["draws"]:
            for shift, width in zip(draw["edge_shifts_eV"], (0.01, 0.02, 0.03)):
                self.assertLessEqual(abs(shift), width + 0.01)

    def test_uniform_common_shift_scales_flux_but_not_populations(self):
        result = propagate_barrier_uncertainty(
            self.rates, (0.0, 0.0, 0.0), samples=10, seed=2,
            shared_half_width_eV=0.04, **self.label
        )
        baseline = result["baseline"]
        for draw in result["draws"]:
            expected = baseline["flux_per_site_s"] * math.exp(-draw["common_shift_eV"] / (KB_EV * 298.15))
            self.assertAlmostEqual(draw["flux_per_site_s"], expected, places=12)
            for p, p0 in zip(draw["populations"], baseline["populations"]):
                self.assertAlmostEqual(p, p0, places=13)

    def test_rate_control_summation_and_step_convergence(self):
        # Scaling every common TS scales all rates and net flux equally:
        # Euler homogeneity independently requires sum of X_i = 1.
        coarse = transition_state_rate_control(self.rates, perturbation_eV=1e-5, **self.label)
        fine = transition_state_rate_control(self.rates, perturbation_eV=5e-6, **self.label)
        self.assertAlmostEqual(fine["coefficient_sum"], 1.0, places=7)
        for x, y in zip(coarse["rate_control_coefficients"], fine["rate_control_coefficients"]):
            self.assertAlmostEqual(x, y, places=7)

    def test_negative_net_flux_control_has_same_consistent_definition(self):
        reversed_driving = (0.01, 0.5, 0.8, 3.0, 0.03, 0.06)
        result = transition_state_rate_control(reversed_driving, **self.label)
        self.assertLess(result["baseline"]["flux_per_site_s"], 0.0)
        self.assertAlmostEqual(result["coefficient_sum"], 1.0, places=7)

    def test_unsupported_physical_label_and_missing_label_are_rejected(self):
        with self.assertRaises(TypeError):
            evaluate_network(self.rates)
        for label in (None, "physical_approved", "LITERATURE FACT", "EXECUTED RESULT", ""):
            with self.assertRaises(ValueError):
                evaluate_network(self.rates, evidence_label=label)
            with self.assertRaises(ValueError):
                transition_state_rate_control(self.rates, evidence_label=label)
            with self.assertRaises(ValueError):
                propagate_barrier_uncertainty(self.rates, (0.01, 0.01, 0.01), evidence_label=label, samples=1, seed=1)

    def test_input_validation_and_sample_limit(self):
        for rates in ((1, 1), (1, 1, 1, 1, 1, 0), (1, 1, 1, 1, 1, math.nan)):
            with self.assertRaises(ValueError):
                evaluate_network(rates, **self.label)
        for sample_count in (0, 10_001, 1.5, True):
            with self.assertRaises(ValueError):
                propagate_barrier_uncertainty(self.rates, (0, 0, 0), samples=sample_count, seed=1, **self.label)
        for widths in ((0.1, -0.1, 0.1), (0.1,), (math.inf, 0.1, 0.1)):
            with self.assertRaises(ValueError):
                propagate_barrier_uncertainty(self.rates, widths, samples=1, seed=1, **self.label)
        for temperature in (0, -1, math.nan, 5e-324):
            with self.assertRaises(ValueError):
                perturb_transition_state(self.rates, 0, 0.1, temperature_K=temperature, **self.label)


if __name__ == "__main__":
    unittest.main()
