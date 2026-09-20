"""Neutral A(sol)+* <=> A* <=> P* <=> P(sol)+*.

Conventions
-----------
All state and transition-state energies are standard Gibbs energies in eV per
molecule on one common reference. Solution activity is c/c_standard (ideal,
dilute-solution assumption); all returned rate constants have units s^-1 because
adsorption multiplies a dimensionless activity. One common Eyring prefactor is
used with transmission coefficient 1. Rates are a verification example, not a
validated reaction model. Site density is mol of sites per geometric m2, and flux
is mol/(geometric m2 s). The single-site model assumes no lateral interactions.

Film transport is an effective stagnant-boundary-layer approximation, NOT a
selective ion membrane or Nernst-Planck model. The plug-flow kernel is isothermal,
constant volumetric flow, one fluid phase, uniform accessible site density and
area, no dispersion, ohmic loss, charging or migration. It cannot represent the
Cu reduction/coupling/ring-opening chemistry without a new validated network.
"""

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

K_B_EV_K = 8.617333262145e-5
K_B_J_K = 1.380649e-23
PLANCK_J_S = 6.62607015e-34
FARADAY_C_MOL = 96485.33212


def _finite(*values: float) -> None:
    if not all(math.isfinite(float(x)) for x in values):
        raise ValueError("All inputs must be finite.")


@dataclass(frozen=True)
class CycleEnergies:
    empty_ev: float
    ads_a_ev: float
    ads_p_ev: float
    solution_a_ev: float
    solution_p_ev: float
    ts_ads_ev: float
    ts_convert_ev: float
    ts_des_ev: float
    temperature_k: float = 298.15


@dataclass(frozen=True)
class RateConstants:
    ads_forward_s: float
    ads_reverse_s: float
    convert_forward_s: float
    convert_reverse_s: float
    des_forward_s: float
    des_reverse_s: float
    equilibrium_p_over_a: float
    log_cycle_closure_residual: float


@dataclass(frozen=True)
class SurfaceState:
    theta_empty: float
    theta_a: float
    theta_p: float
    net_tof_s: float
    edge_fluxes_s: tuple[float, float, float]
    site_balance_residual: float
    steady_residual_s: float


@dataclass(frozen=True)
class FilmState:
    flux_mol_m2_s: float
    surface_a_mol_m3: float
    surface_p_mol_m3: float
    surface: SurfaceState
    flux_balance_residual_mol_m2_s: float


@dataclass(frozen=True)
class PlugFlowState:
    axial_fraction: np.ndarray
    a_mol_m3: np.ndarray
    p_mol_m3: np.ndarray
    product_outlet_mol_s: float
    net_product_formation_mol_s: float
    conversion: float | None
    maximum_material_balance_residual_mol_m3: float
    solver: str


@dataclass(frozen=True)
class ElectricalAccounting:
    energy_j: float
    energy_kwh: float
    energy_kwh_per_kg_isolated: float
    charge_c: float


@dataclass(frozen=True)
class InputEvidence:
    evidence_kind: str
    source_or_run_id: str
    domain: str
    validated: bool
    parameter_use: str


def rate_constants(energies: CycleEnergies) -> RateConstants:
    """Derive all six rates from the SAME states and three transition states.

Negative standard-state activation free energies and numerical under/overflow
are rejected. Such cases require a different physical rate law, not clipping.
"""
    e = energies
    _finite(*e.__dict__.values())
    if e.temperature_k <= 0:
        raise ValueError("Temperature must be positive.")
    endpoints = (
        (e.empty_ev + e.solution_a_ev, e.ads_a_ev, e.ts_ads_ev),
        (e.ads_a_ev, e.ads_p_ev, e.ts_convert_ev),
        (e.ads_p_ev, e.empty_ev + e.solution_p_ev, e.ts_des_ev),
    )
    log_prefactor = math.log(K_B_J_K * e.temperature_k / PLANCK_J_S)
    kt = K_B_EV_K * e.temperature_k
    logs = []
    for initial, final, ts in endpoints:
        if ts < max(initial, final):
            raise ValueError("Transition state lies below an endpoint.")
        logs.extend([log_prefactor - (ts - initial) / kt,
                     log_prefactor - (ts - final) / kt])
    log_equilibrium = -(e.solution_p_ev - e.solution_a_ev) / kt
    if min(logs + [log_equilibrium]) < -700 or max(logs + [log_equilibrium]) > 700:
        raise ValueError("Numerical range exceeded; use a log-domain solver.")
    log_cycle_ratio = logs[0] - logs[1] + logs[2] - logs[3] + logs[4] - logs[5]
    return RateConstants(*map(math.exp, logs), math.exp(log_equilibrium),
                         log_cycle_ratio - log_equilibrium)


def steady_surface(rates: RateConstants, activity_a: float, activity_p: float) -> SurfaceState:
    """Solve the linear stationary site probabilities by matrix-tree weights.

The positive tree formula is algebraically independent of a linear-system solve,
which is used as a verification oracle in the tests. Scaling prevents overflow.
"""
    _finite(activity_a, activity_p, *rates.__dict__.values())
    if min(activity_a, activity_p) < 0:
        raise ValueError("Activities cannot be negative.")
    constants = tuple(rates.__dict__.values())[:6]
    if min(constants) <= 0:
        raise ValueError("Reversible rate constants must be strictly positive.")
    if rates.equilibrium_p_over_a <= 0:
        raise ValueError("Equilibrium constant must be strictly positive.")
    logs = [math.log(x) for x in constants]
    closure = logs[0] - logs[1] + logs[2] - logs[3] + logs[4] - logs[5] - math.log(rates.equilibrium_p_over_a)
    if abs(closure) > 1e-10 or abs(rates.log_cycle_closure_residual) > 1e-10:
        raise ValueError("Imported rates violate thermodynamic cycle closure.")
    a = rates.ads_forward_s * activity_a       # empty -> A*
    b = rates.ads_reverse_s                    # A* -> empty
    c = rates.convert_forward_s                # A* -> P*
    d = rates.convert_reverse_s                # P* -> A*
    e = rates.des_forward_s                    # P* -> empty
    f = rates.des_reverse_s * activity_p        # empty -> P*
    _finite(a, b, c, d, e, f)
    scale = max(a, b, c, d, e, f)
    aa, bb, cc, dd, ee, ff = np.array([a, b, c, d, e, f]) / scale
    weights = np.array([bb * ee + bb * dd + cc * ee,
                        aa * dd + aa * ee + ff * dd,
                        ff * cc + ff * bb + aa * cc])
    if not np.isfinite(weights.sum()) or weights.sum() <= 0:
        raise ValueError("Disconnected or numerically degenerate state graph.")
    t0, ta, tp = weights / weights.sum()
    fluxes = (a * t0 - b * ta, c * ta - d * tp, e * tp - f * t0)
    return SurfaceState(float(t0), float(ta), float(tp), float(np.mean(fluxes)),
                        tuple(map(float, fluxes)), float(t0 + ta + tp - 1),
                        float(max(fluxes) - min(fluxes)))


def film_flux(rates: RateConstants, bulk_a_mol_m3: float, bulk_p_mol_m3: float,
              site_density_mol_m2: float, km_a_m_s: float, km_p_m_s: float,
              standard_concentration_mol_m3: float = 1000.0) -> FilmState:
    """Jointly solve intrinsic kinetics and A/P film transport by one flux root.

J = km_A(c_A,bulk-c_A,surf) = km_P(c_P,surf-c_P,bulk)
  = site_density * TOF(c_A,surf/c_standard, c_P,surf/c_standard).
Reverse reaction gives negative J and is allowed. No electron flux is inferred.
"""
    ca, cp, gamma, ka, kp, cstd = (bulk_a_mol_m3, bulk_p_mol_m3,
        site_density_mol_m2, km_a_m_s, km_p_m_s, standard_concentration_mol_m3)
    _finite(ca, cp, gamma, ka, kp, cstd)
    if min(ca, cp, gamma) < 0 or min(ka, kp, cstd) <= 0:
        raise ValueError("Nonnegative concentrations/sites and positive km/cstd required.")

    def evaluate(flux):
        # Tiny round-off at the physical bracket endpoints is clipped, no rate is.
        sa, sp = max(0.0, ca - flux / ka), max(0.0, cp + flux / kp)
        state = steady_surface(rates, sa / cstd, sp / cstd)
        return sa, sp, state

    lower, upper = -kp * cp, ka * ca
    _finite(lower, upper)
    flux_scale = max(abs(lower), abs(upper))
    if flux_scale == 0 or gamma == 0:
        flux = 0.0
    else:
        def residual(scaled_flux):
            flux = scaled_flux * flux_scale
            result = (flux - gamma * evaluate(flux)[2].net_tof_s) / flux_scale
            _finite(result)
            return result
        flux = flux_scale * brentq(residual, lower / flux_scale, upper / flux_scale,
                                   xtol=1e-14, rtol=1e-13)
    sa, sp, state = evaluate(flux)
    return FilmState(flux, sa, sp, state, flux - gamma * state.net_tof_s)


def _flow_inputs(inlet_a, inlet_p, volume, flow, area_density, cells):
    _finite(inlet_a, inlet_p, volume, flow, area_density)
    if min(inlet_a, inlet_p, volume, area_density) < 0 or flow <= 0:
        raise ValueError("Nonnegative inlet/volume/area and positive flow required.")
    if not isinstance(cells, int) or cells < 1:
        raise ValueError("cells must be a positive integer.")


def _flow_result(a, inlet_a, inlet_p, flow, solver):
    total = inlet_a + inlet_p
    p = total - a
    return PlugFlowState(np.linspace(0, 1, len(a)), a, p, float(flow * p[-1]),
                         float(flow * (p[-1] - inlet_p)),
                         None if inlet_a == 0 else float(1 - a[-1] / inlet_a),
                         float(np.max(np.abs(a + p - total))), solver)


def plug_flow(rates: RateConstants, inlet_a_mol_m3: float, inlet_p_mol_m3: float,
              reactor_volume_m3: float, flow_m3_s: float,
              area_density_m2_m3: float, site_density_mol_m2: float,
              km_a_m_s: float, km_p_m_s: float, cells: int = 100,
              standard_concentration_mol_m3: float = 1000.0) -> PlugFlowState:
    """Positive conservative backward-Euler finite-volume steady plug flow.

Geometric reactive area = reactor volume * area density. A/P inventory is
conserved exactly by the stoichiometric extent coordinate. Grid error is FIRST
ORDER and must be checked independently. Negative conversion denotes net reverse
reaction and is retained. Formed product excludes product fed at the inlet.
"""
    ca, cp = inlet_a_mol_m3, inlet_p_mol_m3
    _flow_inputs(ca, cp, reactor_volume_m3, flow_m3_s, area_density_m2_m3, cells)
    total = ca + cp
    factor = reactor_volume_m3 * area_density_m2_m3 / (flow_m3_s * cells)
    _finite(total, factor)
    profile = [ca]
    # Validate film arguments even for a zero-volume reactor.
    film_flux(rates, ca, cp, site_density_mol_m2, km_a_m_s, km_p_m_s,
              standard_concentration_mol_m3)
    for _ in range(cells):
        prev = profile[-1]
        if total == 0 or factor == 0:
            profile.append(prev)
            continue
        def residual(fraction_a):
            next_a = total * fraction_a
            j = film_flux(rates, next_a, total - next_a, site_density_mol_m2,
                          km_a_m_s, km_p_m_s, standard_concentration_mol_m3).flux_mol_m2_s
            return (next_a - prev + factor * j) / total
        fraction = brentq(residual, 0.0, 1.0, xtol=1e-13, rtol=1e-13)
        profile.append(total * fraction)
    return _flow_result(np.array(profile), ca, cp, flow_m3_s, "backward_euler")


def plug_flow_adaptive(rates: RateConstants, inlet_a_mol_m3: float, inlet_p_mol_m3: float,
                       reactor_volume_m3: float, flow_m3_s: float,
                       area_density_m2_m3: float, site_density_mol_m2: float,
                       km_a_m_s: float, km_p_m_s: float, cells: int = 100,
                       standard_concentration_mol_m3: float = 1000.0) -> PlugFlowState:
    """Independent adaptive Radau axial integrator for grid verification.

This implements the same physical approximations, not independent chemistry.
"""
    ca, cp = inlet_a_mol_m3, inlet_p_mol_m3
    _flow_inputs(ca, cp, reactor_volume_m3, flow_m3_s, area_density_m2_m3, cells)
    total = ca + cp
    factor = reactor_volume_m3 * area_density_m2_m3 / flow_m3_s
    _finite(total, factor)
    film_flux(rates, ca, cp, site_density_mol_m2, km_a_m_s, km_p_m_s,
              standard_concentration_mol_m3)
    if total == 0 or factor == 0:
        return _flow_result(np.full(cells + 1, ca), ca, cp, flow_m3_s, "radau")
    def derivative(_x, y):
        # The invariant interval is enforced only to evaluate intermediate solver
        # stages; the final profile is checked against it below.
        a = float(np.clip(y[0], 0, 1)) * total
        j = film_flux(rates, a, total - a, site_density_mol_m2, km_a_m_s, km_p_m_s,
                      standard_concentration_mol_m3).flux_mol_m2_s
        return [-factor * j / total]
    solution = solve_ivp(derivative, (0, 1), [ca / total], method="Radau",
                         t_eval=np.linspace(0, 1, cells + 1), rtol=1e-10, atol=1e-12)
    if not solution.success:
        raise RuntimeError(solution.message)
    if np.min(solution.y) < -1e-9 or np.max(solution.y) > 1 + 1e-9:
        raise RuntimeError("Adaptive integration left the physical invariant interval.")
    return _flow_result(solution.y[0] * total, ca, cp, flow_m3_s, "radau")


def cell_energy(time_s: Sequence[float], full_cell_voltage_v: Sequence[float],
                total_current_a: Sequence[float], isolated_product_mass_kg: float
                ) -> ElectricalAccounting:
    """Trapezoidal integral of measured/supplied full-cell U*I, never overpotential.

Positive voltage and current mean electrical input. These arrays must come from
the entire cell on one time base. The caller must check sampling convergence;
non-Faradaic current remains included in energy and charge. Auxiliary loads and
separation energy are excluded and must be accounted for separately.
"""
    t, v, i = [np.asarray(x, dtype=float) for x in
               (time_s, full_cell_voltage_v, total_current_a)]
    if any(x.ndim != 1 for x in (t, v, i)) or len(t) < 2 or not (t.shape == v.shape == i.shape):
        raise ValueError("Matching one-dimensional arrays of length >=2 required.")
    if not all(np.all(np.isfinite(x)) for x in (t, v, i)):
        raise ValueError("Nonfinite electrical input.")
    _finite(isolated_product_mass_kg)
    if np.any(np.diff(t) <= 0) or np.any(v < 0) or np.any(i < 0) or isolated_product_mass_kg <= 0:
        raise ValueError("Increasing time, nonnegative input U/I, and positive mass required.")
    # np.trapezoid was introduced in NumPy 2; support NumPy 1.24 as well.
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    try:
        with np.errstate(over="raise", invalid="raise"):
            energy_j = float(integrate(v * i, t))
            charge_c = float(integrate(i, t))
    except FloatingPointError as exc:
        raise ValueError("Electrical input exceeds numerical range.") from exc
    energy_kwh = energy_j / 3_600_000
    _finite(energy_j, charge_c, energy_kwh / isolated_product_mass_kg)
    return ElectricalAccounting(energy_j, energy_kwh,
                                 energy_kwh / isolated_product_mass_kg, charge_c)


def faradaic_efficiency(formed_product_mol: float, charge_c: float,
                        electrons_per_product: int) -> float:
    """FE from formed product BEFORE separation and validated electron stoichiometry.

The neutral verification kernel does not supply electron stoichiometry. Never
derive it from this A->P test model. Values above one are retained only as an
explicit error; investigate analysis, background subtraction and stoichiometry.
"""
    _finite(formed_product_mol, charge_c)
    if formed_product_mol < 0 or charge_c <= 0:
        raise ValueError("Nonnegative formed product and positive charge required.")
    if isinstance(electrons_per_product, bool) or not isinstance(electrons_per_product, int) or electrons_per_product <= 0:
        raise ValueError("Known positive integer electron stoichiometry required.")
    fe = formed_product_mol * electrons_per_product * FARADAY_C_MOL / charge_c
    if fe > 1 + 1e-12:
        raise ValueError("FE exceeds unity; input or stoichiometry is inconsistent.")
    return fe


def target_ranking_gate(evidence: InputEvidence) -> dict:
    """Fail closed: this neutral kernel cannot issue target-electrosynthesis ranks.

Even real neutral reaction parameters cannot certify a Cu electroreduction
prediction. A future target network and audited interface/stoichiometry must have
their own reviewed acceptance gate; changing metadata cannot unlock this kernel.
"""
    reasons = ["neutral_unimolecular_kernel_is_not_target_electrosynthesis_network"]
    if evidence.parameter_use == "unit_fixture":
        reasons.append("synthetic_fixture_cannot_support_physical_prediction")
    if evidence.evidence_kind not in {"LITERATURE", "CALCULATION"}:
        reasons.append("parameter_evidence_is_not_literature_or_executed_calculation")
    if not evidence.source_or_run_id.strip():
        reasons.append("missing_source_or_run_id")
    if not evidence.validated:
        reasons.append("parameters_not_validated")
    if evidence.domain != "target_reaction":
        reasons.append("outside_target_reaction_domain")
    return {"target_ranking_allowed": False, "rejection_reasons": reasons,
            "permitted_use": "software_verification_only",
            "scientific_maturity": "A_engineering_prototype"}
