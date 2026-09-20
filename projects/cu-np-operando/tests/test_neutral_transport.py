"""Software checks using synthetic neutral A/P fixtures, NEVER Cu parameters."""
import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from neutral_transport import (
    CycleEnergies, InputEvidence, rate_constants, steady_surface, film_flux,
    plug_flow, plug_flow_adaptive, cell_energy, faradaic_efficiency,
    target_ranking_gate, K_B_EV_K, FARADAY_C_MOL,
)


def equal_fixture():
    """SYNTHETIC: three equal TS and isoenergetic states, no chemical identity."""
    return rate_constants(CycleEnergies(0, 0, 0, 0, 0, 0.65, 0.65, 0.65))


def flow_inputs():
    return dict(inlet_a_mol_m3=100.0, inlet_p_mol_m3=0.0,
                reactor_volume_m3=1e-5, flow_m3_s=1e-8,
                area_density_m2_m3=2000.0, site_density_mol_m2=1e-5,
                km_a_m_s=1e-6, km_p_m_s=1e-6)


def exact_equal_outlet(rates, inputs):
    """Independent analytic solution, equal kinetic rates and equal film km.

With total concentration Ct fixed, J=lambda*(As-Ps), where lambda =
    Gamma*k/[3(cstd+Ct)]. The coupled bulk flux is
    J=lambda/(1+2*lambda/km)*(Ab-Pb).
The steady PFR consequently has exponential approach to Ct/2.
"""
    total = inputs["inlet_a_mol_m3"] + inputs["inlet_p_mol_m3"]
    intrinsic = inputs["site_density_mol_m2"] * rates.ads_forward_s / (3 * (1000 + total))
    effective = intrinsic / (1 + 2 * intrinsic / inputs["km_a_m_s"])
    area = inputs["area_density_m2_m3"] * inputs["reactor_volume_m3"]
    exponent = -2 * area * effective / inputs["flow_m3_s"]
    return total / 2 + (inputs["inlet_a_mol_m3"] - total / 2) * math.exp(exponent)


class Thermodynamics(unittest.TestCase):
    def test_imported_unphysical_rates_rejected(self):
        from neutral_transport.core import RateConstants
        for invalid in (RateConstants(1, 1, 1, 1, 1, 1, -1, 100),
                        RateConstants(2, 1, 1, 1, 1, 1, 1, 0),
                        RateConstants(1, 1, 1, 1, 1, 1, 1, 100)):
            with self.assertRaises(ValueError):
                steady_surface(invalid, 1.0, 1.0)

    def test_closed_thermodynamic_cycle(self):
        energies = CycleEnergies(0.0, -0.12, -0.08, 0.0, -0.04, 0.65, 0.58, 0.62)
        rates = rate_constants(energies)
        self.assertLess(abs(rates.log_cycle_closure_residual), 1e-12)
        product_ratio = (rates.ads_forward_s / rates.ads_reverse_s
                         * rates.convert_forward_s / rates.convert_reverse_s
                         * rates.des_forward_s / rates.des_reverse_s)
        self.assertAlmostEqual(product_ratio / rates.equilibrium_p_over_a, 1, places=12)

    def test_equilibrium_has_zero_edge_flux(self):
        energies = CycleEnergies(0.0, -0.12, -0.08, 0.0, -0.04, 0.65, 0.58, 0.62)
        rates = rate_constants(energies)
        state = steady_surface(rates, 0.01, 0.01 * rates.equilibrium_p_over_a)
        self.assertLess(max(abs(j) for j in state.edge_fluxes_s), 1e-12)
        expected_a_over_empty = (0.01 * math.exp(
            -(energies.ads_a_ev - energies.empty_ev - energies.solution_a_ev)
            / (K_B_EV_K * energies.temperature_k)))
        self.assertAlmostEqual(state.theta_a / state.theta_empty, expected_a_over_empty, places=12)

    def test_equal_rates_analytic_tof(self):
        rates = equal_fixture()
        for a, p in [(0, 0), (0.1, 0), (0.1, 0.2), (10, 3), (1, 1)]:
            state = steady_surface(rates, a, p)
            exact = rates.ads_forward_s * (a - p) / (3 * (1 + a + p))
            self.assertAlmostEqual(state.net_tof_s, exact, places=11)
            self.assertLess(abs(state.site_balance_residual), 1e-14)
            self.assertLess(state.steady_residual_s, 1e-11)
            self.assertGreaterEqual(min(state.theta_empty, state.theta_a, state.theta_p), 0)

    def test_independent_linear_system(self):
        rates = rate_constants(CycleEnergies(0, -0.03, -0.02, 0, -0.01, 0.60, 0.62, 0.64))
        aa, pp = 0.02, 0.003
        a, b = rates.ads_forward_s * aa, rates.ads_reverse_s
        c, d = rates.convert_forward_s, rates.convert_reverse_s
        e, f = rates.des_forward_s, rates.des_reverse_s * pp
        matrix = np.array([[-a-f, b, e], [a, -b-c, d], [1, 1, 1]])
        direct = np.linalg.solve(matrix, np.array([0, 0, 1]))
        state = steady_surface(rates, aa, pp)
        np.testing.assert_allclose(direct, [state.theta_empty, state.theta_a, state.theta_p], rtol=1e-12)

    def test_invalid_barrier_and_temperature_rejected(self):
        with self.assertRaises(ValueError):
            rate_constants(CycleEnergies(0, 1, 0, 0, 0, 0.5, 1.2, 0.5))
        with self.assertRaises(ValueError):
            rate_constants(CycleEnergies(0, 0, 0, 0, 0, 1, 1, 1, temperature_k=0))
        with self.assertRaises(ValueError):
            steady_surface(equal_fixture(), -1, 0)


class Transport(unittest.TestCase):
    def test_joint_film_solution_matches_analytic(self):
        rates = equal_fixture()
        gamma, km, ca, cp = 1e-5, 1e-6, 100, 3
        state = film_flux(rates, ca, cp, gamma, km, km)
        intrinsic = gamma * rates.ads_forward_s / (3 * (1000 + ca + cp))
        exact = intrinsic * (ca - cp) / (1 + 2 * intrinsic / km)
        self.assertAlmostEqual(state.flux_mol_m2_s / exact, 1, places=11)
        self.assertLess(abs(state.flux_balance_residual_mol_m2_s), 1e-16)
        self.assertAlmostEqual(km * (ca - state.surface_a_mol_m3) / exact, 1, places=11)
        self.assertAlmostEqual(km * (state.surface_p_mol_m3 - cp) / exact, 1, places=11)

    def test_transport_asymptotic_limits(self):
        rates = equal_fixture()
        ca, cp, gamma = 100, 0, 1e-5
        fast = film_flux(rates, ca, cp, gamma, 1e3, 1e3)
        intrinsic = gamma * steady_surface(rates, ca / 1000, 0).net_tof_s
        self.assertAlmostEqual(fast.flux_mol_m2_s / intrinsic, 1, places=8)
        # Very slow equal A/P films: approach two-film equilibrium limit km*CA/2.
        km = 1e-12
        slow = film_flux(rates, ca, cp, gamma, km, km)
        self.assertAlmostEqual(slow.flux_mol_m2_s / (km * ca / 2), 1, places=5)

    def test_reverse_reaction_and_zero_sites(self):
        reverse = film_flux(equal_fixture(), 1, 100, 1e-5, 1e-6, 2e-6)
        self.assertLess(reverse.flux_mol_m2_s, 0)
        self.assertGreaterEqual(reverse.surface_a_mol_m3, 0)
        self.assertGreaterEqual(reverse.surface_p_mol_m3, 0)
        zero = film_flux(equal_fixture(), 100, 0, 0, 1e-6, 1e-6)
        self.assertEqual(zero.flux_mol_m2_s, 0)
        self.assertEqual(zero.surface_a_mol_m3, 100)

    def test_no_reaction_or_zero_volume(self):
        for delta in [{"site_density_mol_m2": 0}, {"reactor_volume_m3": 0},
                      {"inlet_a_mol_m3": 0}]:
            params = flow_inputs() | delta
            result = plug_flow(equal_fixture(), **params, cells=10)
            np.testing.assert_array_equal(result.a_mol_m3, params["inlet_a_mol_m3"])
            self.assertEqual(result.net_product_formation_mol_s, 0)

    def test_pfr_grid_converges_to_independent_analytic_limit(self):
        rates, params = equal_fixture(), flow_inputs()
        exact = exact_equal_outlet(rates, params)
        errors = []
        for grid in [20, 40, 80, 160]:
            result = plug_flow(rates, **params, cells=grid)
            errors.append(abs(result.a_mol_m3[-1] - exact))
            self.assertLess(result.maximum_material_balance_residual_mol_m3, 1e-12)
            self.assertGreaterEqual(float(result.a_mol_m3.min()), 0)
            self.assertGreaterEqual(float(result.p_mol_m3.min()), 0)
        ratios = np.array(errors[:-1]) / errors[1:]
        self.assertTrue(np.all((ratios > 1.8) & (ratios < 2.2)), ratios)
        self.assertLess(errors[-1] / params["inlet_a_mol_m3"], 0.001)

    def test_adaptive_independent_integrator(self):
        rates, params = equal_fixture(), flow_inputs()
        exact = exact_equal_outlet(rates, params)
        adaptive = plug_flow_adaptive(rates, **params)
        self.assertLess(abs(adaptive.a_mol_m3[-1] - exact), 1e-7)
        grid = plug_flow(rates, **params, cells=160)
        self.assertLess(abs(grid.a_mol_m3[-1] - adaptive.a_mol_m3[-1]) / 100, 0.001)

    def test_nonzero_product_inlet_not_counted_as_formation(self):
        params = flow_inputs() | {"inlet_p_mol_m3": 5}
        result = plug_flow(equal_fixture(), **params)
        inlet_product = 5 * params["flow_m3_s"]
        self.assertAlmostEqual(result.product_outlet_mol_s - result.net_product_formation_mol_s,
                               inlet_product, places=15)
        self.assertAlmostEqual(result.net_product_formation_mol_s,
            params["flow_m3_s"] * (100 - result.a_mol_m3[-1]), places=15)

    def test_integrated_surface_source_matches_outlet_formation_in_both_directions(self):
        rates = equal_fixture()
        for ca, cp in [(100.0, 5.0), (1.0, 100.0)]:
            params = flow_inputs() | {"inlet_a_mol_m3": ca, "inlet_p_mol_m3": cp}
            cells = 80
            result = plug_flow(rates, **params, cells=cells)
            # Integrate independently evaluated local source using the same
            # right-face quadrature as the backward-Euler discretization.
            source = np.array([film_flux(rates, a, p, params["site_density_mol_m2"],
                        params["km_a_m_s"], params["km_p_m_s"]).flux_mol_m2_s
                       for a, p in zip(result.a_mol_m3[1:], result.p_mol_m3[1:])])
            integrated = (params["reactor_volume_m3"] * params["area_density_m2_m3"]
                          * float(source.sum()) / cells)
            self.assertAlmostEqual(integrated / result.net_product_formation_mol_s, 1, places=9)
            self.assertEqual(np.sign(result.net_product_formation_mol_s), np.sign(ca - cp))


class AccountingAndEvidence(unittest.TestCase):
    def test_full_cell_energy_and_si_units(self):
        result = cell_energy([0, 1800, 3600], [3, 3, 3], [2, 2, 2], 0.001)
        self.assertAlmostEqual(result.energy_j, 21600)
        self.assertAlmostEqual(result.energy_kwh, 0.006)
        self.assertAlmostEqual(result.energy_kwh_per_kg_isolated, 6)
        self.assertAlmostEqual(result.charge_c, 7200)

    def test_variable_power_integration(self):
        # Constant I and linearly increasing U: trapezoid is analytically exact.
        result = cell_energy([0, 2, 5], [1, 3, 6], [2, 2, 2], 0.01)
        self.assertAlmostEqual(result.energy_j, 35)
        self.assertAlmostEqual(result.charge_c, 10)

    def test_fe_is_formed_moles_not_isolated_mass(self):
        fe = faradaic_efficiency(0.01, 7200, 2)
        self.assertAlmostEqual(fe, 0.02 * FARADAY_C_MOL / 7200)
        high_recovery = cell_energy([0, 3600], [3, 3], [2, 2], 0.002)
        low_recovery = cell_energy([0, 3600], [3, 3], [2, 2], 0.001)
        self.assertAlmostEqual(low_recovery.energy_kwh_per_kg_isolated,
                               2 * high_recovery.energy_kwh_per_kg_isolated)
        self.assertAlmostEqual(fe, faradaic_efficiency(0.01, high_recovery.charge_c, 2))

    def test_invalid_electrical_inputs_are_rejected(self):
        for args in [([0, 0], [3, 3], [2, 2], 1), ([0, 1], [3, 3], [-2, 2], 1),
                     ([0, 1], [3], [2, 2], 1), ([0, 1], [3, 3], [2, 2], 0)]:
            with self.assertRaises(ValueError):
                cell_energy(*args)
        for args in [(1, 1, 2), (0.01, 0, 2), (0.01, 10000, 0)]:
            with self.assertRaises(ValueError):
                faradaic_efficiency(*args)

    def test_fixture_cannot_be_mislabeled_as_physical_rank(self):
        evidence = InputEvidence("HYPOTHESIS", "synthetic-neutral-fixture", "neutral_fixture", False, "unit_fixture")
        result = target_ranking_gate(evidence)
        self.assertFalse(result["target_ranking_allowed"])
        self.assertIn("synthetic_fixture_cannot_support_physical_prediction", result["rejection_reasons"])
        # Metadata relabeling cannot unlock a physically inappropriate network.
        forged = InputEvidence("CALCULATION", "apparently-real-run", "target_reaction", True, "prediction")
        self.assertFalse(target_ranking_gate(forged)["target_ranking_allowed"])

    def test_extreme_finite_inputs_fail_without_nan_outputs(self):
        with self.assertRaises(ValueError):
            steady_surface(equal_fixture(), 1e308, 0)
        with self.assertRaises(ValueError):
            film_flux(equal_fixture(), 1e308, 0, 1e-5, 1e308, 1e-6)
        with self.assertRaises(ValueError):
            plug_flow(equal_fixture(), **(flow_inputs() | {"reactor_volume_m3": 1e308}))
        with self.assertRaises(ValueError):
            cell_energy([0, 1], [1e308, 1e308], [1e308, 1e308], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
