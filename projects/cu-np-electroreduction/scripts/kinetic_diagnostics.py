"""SYNTHETIC SOFTWARE TESTS ONLY: three-state kinetic diagnostics.

This module does not describe the Cu/quinazolinone/nitrile chemistry. Rates
are six effective first-order toy rates ordered 0->1, 1->0, 1->2, 2->1,
2->0, 0->2. A common transition-state perturbation multiplies BOTH rates
of an edge. It preserves that edge's equilibrium ratio and cycle affinity.
An initially driven toy network stays driven; this does not turn its
non-equilibrium steady state into thermodynamic equilibrium.

Only the explicit evidence_label='synthetic_software_test' is accepted.
There is deliberately no physical_approved interface or override.
"""

import math
import random

try:
    from .science import KB_EV, site_steady_state
except ImportError:  # Direct import when scripts/ is placed on sys.path.
    from science import KB_EV, site_steady_state


SYNTHETIC_LABEL = "synthetic_software_test"
SCOPE = "three_state_software_validation_not_target_reaction"
MAX_SAMPLES = 10_000


def _label(value):
    if value != SYNTHETIC_LABEL:
        raise ValueError(
            "Only explicit synthetic_software_test inputs are allowed; "
            "physical_approved is unsupported and requires a separate audited model."
        )


def _finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + " must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(name + " must be finite")
    return value


def _rates(values):
    try:
        values = tuple(values)
    except TypeError as exc:
        raise ValueError("Six finite positive rates are required") from exc
    if len(values) != 6:
        raise ValueError("Six finite positive rates are required")
    result = tuple(_finite_number(x, "rate") for x in values)
    if min(result) <= 0.0:
        raise ValueError("Rates must be positive")
    return result


def _thermal_energy(temperature_K):
    temperature = _finite_number(temperature_K, "temperature_K")
    if temperature <= 0.0:
        raise ValueError("temperature_K must be positive")
    kt = KB_EV * temperature
    if not math.isfinite(kt) or kt <= 0.0:
        raise ValueError("temperature_K is outside the supported numerical range")
    return kt


def _evaluate(rates):
    # Common rescaling avoids overflow in the tree-weight products used by
    # science.site_steady_state. Fail explicitly outside its numerical domain.
    scale = max(rates)
    scaled = tuple(x / scale for x in rates)
    if min(scaled) == 0.0:
        raise ValueError("Rate dynamic range exceeds supported floating-point range")
    theta, normalized_flux = site_steady_state(scaled)
    if not all(math.isfinite(x) for x in theta):
        raise ValueError("Rate dynamic range gives unresolved site populations")
    edge_fluxes = (
        rates[0] * theta[0] - rates[1] * theta[1],
        rates[2] * theta[1] - rates[3] * theta[2],
        rates[4] * theta[2] - rates[5] * theta[0],
    )
    traffic = max(
        scaled[0] * theta[0] + scaled[1] * theta[1],
        scaled[2] * theta[1] + scaled[3] * theta[2],
        scaled[4] * theta[2] + scaled[5] * theta[0],
    )
    return {
        "evidence_label": SYNTHETIC_LABEL,
        "scope": SCOPE,
        "populations": list(theta),
        "flux_per_site_s": normalized_flux * scale,
        "edge_fluxes_per_site_s": list(edge_fluxes),
        "log_forward_reverse_ratios": [
            math.log(rates[i]) - math.log(rates[i + 1]) for i in (0, 2, 4)
        ],
        "cycle_affinity_over_kT": math.fsum(
            math.log(rates[i]) - math.log(rates[i + 1]) for i in (0, 2, 4)
        ),
        "relative_net_flux": abs(normalized_flux) / traffic if traffic else 0.0,
    }


def evaluate_network(rates, *, evidence_label):
    """Return normalized populations and a toy per-site flux, never current."""
    _label(evidence_label)
    return _evaluate(_rates(rates))


def _shift(rates, shifts_eV, kt):
    shifted = []
    for edge, delta in enumerate(shifts_eV):
        exponent = -delta / kt
        if abs(exponent) > 600.0:
            raise ValueError("TS perturbation exceeds supported numerical range")
        factor = math.exp(exponent)
        shifted.extend((rates[2 * edge] * factor, rates[2 * edge + 1] * factor))
    return _rates(shifted)


def perturb_transition_state(
    rates, edge_index, delta_g_eV, *, evidence_label, temperature_K=298.15
):
    """Raise/lower one common toy TS, keeping state differences unchanged.

    The rate-only toy has no absolute TS energies; shifts are relative
    software-test parameters, not validated physical activation barriers.
    """
    _label(evidence_label)
    if isinstance(edge_index, bool) or not isinstance(edge_index, int) or edge_index not in (0, 1, 2):
        raise ValueError("edge_index must be 0, 1 or 2")
    shifts = [0.0, 0.0, 0.0]
    shifts[edge_index] = _finite_number(delta_g_eV, "delta_g_eV")
    return _shift(_rates(rates), shifts, _thermal_energy(temperature_K))


def transition_state_rate_control(
    rates, *, evidence_label, temperature_K=298.15, perturbation_eV=1e-5,
    relative_flux_tolerance=1e-10
):
    """Central differences: X_i = -kT * d ln|J| / d G_TS_i.

    Both forward and reverse barriers on edge i are perturbed together.
    State energies, equilibrium ratios and any toy cycle driving stay fixed.
    Logarithmic net-flux control is undefined at equilibrium and is rejected.
    Repeating with a smaller step tests finite-difference convergence.
    """
    _label(evidence_label)
    baseline_rates = _rates(rates)
    kt = _thermal_energy(temperature_K)
    step = _finite_number(perturbation_eV, "perturbation_eV")
    tolerance = _finite_number(relative_flux_tolerance, "relative_flux_tolerance")
    if step <= 0.0 or not 0.0 < tolerance < 1.0:
        raise ValueError("Positive step and relative flux tolerance in (0,1) required")
    baseline = _evaluate(baseline_rates)
    if baseline["relative_net_flux"] <= tolerance:
        raise ValueError("Log net-flux rate control is undefined or unresolved near equilibrium")
    coefficients = []
    for edge in range(3):
        positive = [step if i == edge else 0.0 for i in range(3)]
        negative = [-x for x in positive]
        plus = _evaluate(_shift(baseline_rates, positive, kt))
        minus = _evaluate(_shift(baseline_rates, negative, kt))
        for perturbed in (plus, minus):
            if perturbed["relative_net_flux"] <= tolerance:
                raise ValueError("Perturbed net flux is numerically unresolved")
            if math.copysign(1.0, perturbed["flux_per_site_s"]) != math.copysign(1.0, baseline["flux_per_site_s"]):
                raise ValueError("Unexpected flux sign change under common-TS perturbation")
        coefficients.append(-kt * (
            math.log(abs(plus["flux_per_site_s"]))
            - math.log(abs(minus["flux_per_site_s"]))
        ) / (2.0 * step))
    return {
        "evidence_label": SYNTHETIC_LABEL,
        "scope": SCOPE,
        "baseline": baseline,
        "temperature_K": float(temperature_K),
        "perturbation_eV": step,
        "rate_control_coefficients": coefficients,
        "coefficient_sum": math.fsum(coefficients),
        "definition": "-kT d ln(abs(net_flux)) / d common_TS_energy",
    }


def _quantiles(values):
    ordered = sorted(values)
    result = {}
    for name, probability in (("q025", 0.025), ("median", 0.5), ("q975", 0.975)):
        position = probability * (len(ordered) - 1)
        low = int(math.floor(position))
        high = int(math.ceil(position))
        result[name] = ordered[low] + (ordered[high] - ordered[low]) * (position - low)
    return result


def propagate_barrier_uncertainty(
    rates, barrier_half_widths_eV, *, evidence_label, samples, seed,
    temperature_K=298.15, shared_half_width_eV=0.0
):
    """Bounded toy Monte Carlo preserving every forward/reverse rate ratio.

    Each draw uses one uniform common offset plus independent uniform edge
    offsets. The common term represents correlated software-test uncertainty.
    Widths are assumptions, not fitted uncertainties or DFT error estimates.
    Samples are retained to make seed, site balance and ratios inspectable.
    Reported quantiles are conditional toy-distribution quantiles, not
    confidence/credible intervals for the real electrocatalytic reaction.
    """
    _label(evidence_label)
    baseline_rates = _rates(rates)
    kt = _thermal_energy(temperature_K)
    if isinstance(samples, bool) or not isinstance(samples, int) or not 1 <= samples <= MAX_SAMPLES:
        raise ValueError("samples must be an integer between 1 and 10000")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("An explicit integer seed is required")
    try:
        widths = tuple(barrier_half_widths_eV)
    except TypeError as exc:
        raise ValueError("Three nonnegative barrier half-widths are required") from exc
    if len(widths) != 3:
        raise ValueError("Three nonnegative barrier half-widths are required")
    widths = tuple(_finite_number(x, "barrier_half_width") for x in widths)
    shared = _finite_number(shared_half_width_eV, "shared_half_width_eV")
    if min(widths) < 0.0 or shared < 0.0:
        raise ValueError("Uncertainty half-widths cannot be negative")
    if (max(widths) + shared) / kt > 600.0:
        raise ValueError("Uncertainty support exceeds supported numerical range")
    rng = random.Random(seed)
    records = []
    for index in range(samples):
        common_shift = rng.uniform(-shared, shared)
        shifts = [common_shift + rng.uniform(-width, width) for width in widths]
        sampled_rates = _shift(baseline_rates, shifts, kt)
        diagnostic = _evaluate(sampled_rates)
        records.append({
            "sample_index": index,
            "evidence_label": SYNTHETIC_LABEL,
            "common_shift_eV": common_shift,
            "edge_shifts_eV": shifts,
            "rates_per_s": list(sampled_rates),
            **diagnostic,
        })
    return {
        "evidence_label": SYNTHETIC_LABEL,
        "scope": SCOPE,
        "seed": seed,
        "samples": samples,
        "temperature_K": float(temperature_K),
        "barrier_half_widths_eV": list(widths),
        "shared_half_width_eV": shared,
        "distribution": "bounded uniform shared offset plus independent bounded uniform edge offsets",
        "interpretation": "synthetic assumed-distribution quantiles, not target-reaction predictive intervals",
        "baseline": _evaluate(baseline_rates),
        "flux_quantiles_per_site_s": _quantiles([r["flux_per_site_s"] for r in records]),
        "population_quantiles": [_quantiles([r["populations"][i] for r in records]) for i in range(3)],
        "draws": records,
    }
