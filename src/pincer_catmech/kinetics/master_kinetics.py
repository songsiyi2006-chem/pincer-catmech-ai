"""Eighteen-state, eighteen-channel batch network with exact derivatives.

This is a declared mechanistic *hypothesis*, not eighteen established elementary
steps. In particular the bound-alcohol poisoning and arm-cleavage channels need
their own structures and barriers. No rate is inferred from a reaction energy.

Concentrations are M, time is s, standard concentration is 1 M. Free tert-butoxide
is a dynamic species; it is not silently equated with total added potassium salt.
tBuOH is an explicitly declared neutral reservoir through its dimensionless
activity. The irreversible cleavage fragment is eliminated exactly: its amount
formed equals the increase in ``cat_degraded``. This closes the elemental and
charge ledger without a nineteenth ODE. K+ is a conserved spectator.

The hydrogen-containing dimer is formed from TWO hydrogenated catalysts, not
from two hydrogen-free active catalysts. It carries two metal atoms and four
additional hydrogen atoms (two hydrides plus two ligand protons).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Mapping
import hashlib
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import eye, csr_matrix, kron

from .microkinetics import ComputedFreeEnergy, EvidenceError, KB_OVER_H, R_KCAL

SPECIES = (
    "cat", "cat_alcohol", "cat_h2_aldehyde", "cat_h2", "cat_h2_imine",
    "cat_product", "dimer", "cat_co", "cat_degraded", "alcohol",
    "aldehyde", "amine", "hemiaminal", "imine", "product", "water",
    "benzene", "free_base",
)
INDEX = {name: i for i, name in enumerate(SPECIES)}
STEP_NAMES = (
    "alcohol_coordination", "alcohol_h_transfer", "aldehyde_release",
    "imine_coordination", "imine_hydrogenation", "product_release",
    "hemiaminal_addition", "direct_dehydration", "base_assisted_addition",
    "base_assisted_dehydration", "one_alcohol_shuttle", "two_alcohol_shuttle",
    "three_alcohol_shuttle", "hydrogenated_dimerization", "decarbonylation",
    "bound_alcohol_decarbonylation", "arm_cleavage", "bound_alcohol_arm_cleavage",
)
# Reactants/products include kinetic catalysts even when they cancel in S.
REACTIONS = (
    ({"cat": 1, "alcohol": 1}, {"cat_alcohol": 1}),
    ({"cat_alcohol": 1}, {"cat_h2_aldehyde": 1}),
    ({"cat_h2_aldehyde": 1}, {"cat_h2": 1, "aldehyde": 1}),
    ({"cat_h2": 1, "imine": 1}, {"cat_h2_imine": 1}),
    ({"cat_h2_imine": 1}, {"cat_product": 1}),
    ({"cat_product": 1}, {"cat": 1, "product": 1}),
    ({"aldehyde": 1, "amine": 1}, {"hemiaminal": 1}),
    ({"hemiaminal": 1}, {"imine": 1, "water": 1}),
    ({"aldehyde": 1, "amine": 1, "free_base": 1}, {"hemiaminal": 1, "free_base": 1}),
    ({"hemiaminal": 1, "free_base": 1}, {"imine": 1, "water": 1, "free_base": 1}),
    ({"hemiaminal": 1, "tbuoh": 1}, {"imine": 1, "water": 1, "tbuoh": 1}),
    ({"hemiaminal": 1, "tbuoh": 2}, {"imine": 1, "water": 1, "tbuoh": 2}),
    ({"hemiaminal": 1, "tbuoh": 3}, {"imine": 1, "water": 1, "tbuoh": 3}),
    ({"cat_h2": 2}, {"dimer": 1}),
    ({"cat": 1, "aldehyde": 1}, {"cat_co": 1, "benzene": 1}),
    ({"cat_alcohol": 1, "aldehyde": 1}, {"cat_co": 1, "benzene": 1, "alcohol": 1}),
    ({"cat": 1, "free_base": 1}, {"cat_degraded": 1, "cleavage_fragment": 1}),
    ({"cat_alcohol": 1, "free_base": 1}, {"cat_degraded": 1, "cleavage_fragment": 1, "alcohol": 1}),
)
IRREVERSIBLE = (14, 15, 16, 17)
S = np.zeros((18, 18), dtype=float)
for _j, (_left, _right) in enumerate(REACTIONS):
    for _side, _sign in ((_left, -1), (_right, 1)):
        for _name, _count in _side.items():
            if _name in INDEX:
                S[INDEX[_name], _j] += _sign * _count
METAL = np.array([1, 1, 1, 1, 1, 1, 2, 1, 1] + [0] * 9, dtype=float)
BASE_EQUIVALENTS = np.array([0] * 8 + [1] + [0] * 8 + [1], dtype=float)
S.setflags(write=False)
METAL.setflags(write=False)
BASE_EQUIVALENTS.setflags(write=False)
assert np.array_equal(METAL @ S, np.zeros(18))
assert np.array_equal(BASE_EQUIVALENTS @ S, np.zeros(18))


def _finite_vector(value, size, name, *, nonnegative=False):
    a = np.asarray(value, dtype=float)
    if a.shape != (size,) or not np.isfinite(a).all():
        raise ValueError(f"{name} must have {size} finite entries")
    if nonnegative and np.any(a < 0):
        raise ValueError(f"{name} cannot contain negative entries")
    return a


@dataclass(frozen=True)
class MasterThermochemistry:
    """One temperature/protocol, independently certified states and all 18 TSs.

    ``cleavage_fragment`` must identify the actual P-arm fragment/product. For
    the proposed nucleophilic cleavage it is R-OtBu and the degraded catalyst
    is anionic. Aryl and alkyl cleavage cannot share an assumed barrier.
    """
    temperature: float
    states: Mapping[str, ComputedFreeEnergy]
    transition_states: Mapping[str, ComputedFreeEnergy]
    irreversible_sink_approximation: str

    def __post_init__(self):
        # A frozen dataclass alone does not protect dictionaries supplied by a
        # caller. Keep source records stable after rate validation.
        object.__setattr__(self, "states", MappingProxyType(dict(self.states)))
        object.__setattr__(self, "transition_states", MappingProxyType(dict(self.transition_states)))

    def validate(self):
        needed = set(SPECIES) | {"tbuoh", "cleavage_fragment"}
        if set(self.states) != needed or set(self.transition_states) != set(STEP_NAMES):
            raise EvidenceError("All 20 state/reference energies and all 18 matched TS certificates are required")
        if not isinstance(self.irreversible_sink_approximation, str) or not self.irreversible_sink_approximation.strip():
            raise EvidenceError("Declare the irreversible poisoning/cleavage sink approximation")
        for state in self.states.values():
            state.validate(self.temperature)
        for name, ts in self.transition_states.items():
            ts.validate(self.temperature, transition_state=True,
                        dehydrogenation=name == "alcohol_h_transfer")
        records = list(self.states.values()) + list(self.transition_states.values())
        if len({record.protocol for record in records}) != 1:
            raise EvidenceError("All states and TSs require one method, solvent and standard-state protocol")

        def composition(side):
            counts = Counter()
            charge = 0
            for name, coefficient in side.items():
                record = self.states[name]
                for symbol, number in record.composition:
                    counts[symbol] += number * coefficient
                charge += record.charge * coefficient
            return dict(counts), charge

        for name, (left, right) in zip(STEP_NAMES, REACTIONS):
            if composition(left) != composition(right):
                raise EvidenceError(f"Unbalanced elements/charge in {name}")
            ts = self.transition_states[name]
            if composition(left) != (dict(ts.composition), ts.charge):
                raise EvidenceError(f"TS composition/charge mismatch in {name}")

        # The reservoir and cleavage ledger have fixed identities in this model.
        # Atom balance alone would also accept a mislabeled species network.
        for name, expected, charge in (
            ("tbuoh", {"C": 4, "H": 10, "O": 1}, 0),
            ("free_base", {"C": 4, "H": 9, "O": 1}, -1),
        ):
            record = self.states[name]
            if dict(record.composition) != expected or record.charge != charge:
                raise EvidenceError(f"{name} must have its declared tert-butanol/tert-butoxide identity")
        if self.states["cleavage_fragment"].charge != 0:
            raise EvidenceError("The declared R-OtBu cleavage fragment must be neutral")
        metals = {"Mn", "Fe", "Co", "Ru"}
        active = dict(self.states["cat"].composition)
        present = metals.intersection(active)
        if len(present) != 1 or active[next(iter(present))] != 1:
            raise EvidenceError("Active catalyst must contain exactly one supported metal atom")
        metal = next(iter(present))
        for i, name in enumerate(SPECIES):
            counts = dict(self.states[name].composition)
            if counts.get(metal, 0) != METAL[i] or (metals.intersection(counts)-{metal}):
                raise EvidenceError(f"Metal identity/inventory mismatch in {name}")
        if any(metals.intersection(dict(self.states[name].composition)) for name in ("tbuoh", "cleavage_fragment")):
            raise EvidenceError("Reservoir and cleavage fragment cannot contain catalyst metal")

    def _eyring_arrays(self):
        """Numeric 1 M coefficients derived only from validated source energies."""
        self.validate()
        forward, reverse = [], []
        prefactor = KB_OVER_H * self.temperature
        for j, (name, (left, right)) in enumerate(zip(STEP_NAMES, REACTIONS)):
            ts_g = self.transition_states[name].gibbs_kcal_mol
            values = []
            for side in (left, right):
                barrier = ts_g - sum(n * self.states[x].gibbs_kcal_mol for x, n in side.items())
                if not math.isfinite(barrier) or barrier < 0:
                    raise EvidenceError(f"Nonfinite or negative standard-state barrier in {name}; review association/diffusion treatment")
                # C_standard=1 M: numeric factor is one, dimensions remain
                # M**(1 - dynamic_molecularity) / s after absorbing the powers
                # of C_standard belonging to dimensionless tBuOH activity.
                values.append(prefactor * math.exp(-barrier / (R_KCAL * self.temperature)))
            forward.append(values[0])
            reverse.append(0.0 if j in IRREVERSIBLE else values[1])
        return np.array(forward), np.array(reverse)

    def rates(self):
        forward, reverse = self._eyring_arrays()
        return MasterRates(forward, reverse, self.temperature,
                           "computed", self)

    def elemental_charge_ledger(self):
        """Rows count atoms/charge, with the eliminated fragment assigned to X.

        The actual inventory is L @ c minus the constant fragment composition
        times c_X(0), plus any independently specified initial free fragment.
        This constant offset does not affect the conservation residual.
        """
        elements = sorted({symbol for record in self.states.values() for symbol, _ in record.composition})
        labels = tuple(elements) + ("charge",)
        fragment = self.states["cleavage_fragment"]
        matrix = np.array([[dict(self.states[name].composition).get(symbol, 0) for name in SPECIES]
                           for symbol in elements] + [[self.states[name].charge for name in SPECIES]], dtype=float)
        matrix[:, INDEX["cat_degraded"]] += np.array(
            [dict(fragment.composition).get(symbol, 0) for symbol in elements] + [fragment.charge])
        if np.any(matrix @ S != 0):
            raise EvidenceError("Reduced atom/charge ledger is not conserved")
        return labels, matrix


@dataclass(frozen=True)
class MasterRates:
    forward: np.ndarray
    reverse: np.ndarray
    temperature: float
    evidence_kind: str
    thermochemistry: MasterThermochemistry | None = None
    ts_perturbations_theta: tuple[float, ...] = (0.0,) * 18

    def __post_init__(self):
        for name in ("forward", "reverse"):
            array = _finite_vector(getattr(self, name), 18, name, nonnegative=True).copy()
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        if not math.isfinite(self.temperature) or self.temperature <= 0:
            raise ValueError("Positive finite temperature required")
        if np.any(self.reverse[list(IRREVERSIBLE)] != 0):
            raise ValueError("Declared irreversible sinks require zero reverse rates")
        theta = _finite_vector(self.ts_perturbations_theta, 18, "TS perturbations")
        object.__setattr__(self, "ts_perturbations_theta", tuple(float(x) for x in theta))
        if self.evidence_kind not in {"computed", "computed_TS_perturbation", "software_fixture"}:
            raise EvidenceError("Rate evidence must be computed or an explicit software fixture")
        if self.evidence_kind in {"computed", "computed_TS_perturbation"}:
            if self.thermochemistry is None:
                raise EvidenceError("Computed rates require complete thermochemistry")
            if not math.isclose(self.temperature, self.thermochemistry.temperature, rel_tol=0, abs_tol=1e-6):
                raise EvidenceError("Rate and thermochemistry temperatures differ")
            if (self.evidence_kind == "computed") != (not np.any(theta)):
                raise EvidenceError("Nonzero TS energy perturbations require their explicit evidence label")
            forward, reverse = self.thermochemistry._eyring_arrays()
            with np.errstate(over="raise", invalid="raise"):
                try:
                    expected = (forward*np.exp(theta), reverse*np.exp(theta))
                except FloatingPointError as error:
                    raise EvidenceError("TS perturbation exceeds the numeric range") from error
            if any(not np.allclose(actual, reference, rtol=1e-12, atol=0)
                   for actual, reference in zip((self.forward, self.reverse), expected)):
                raise EvidenceError("Computed rate arrays differ from their source Eyring rates/declared TS perturbations")
        elif self.thermochemistry is not None:
            raise EvidenceError("Software fixtures cannot carry computed thermochemistry")

    def perturbed(self, step: int, theta: float):
        """Lower one TS by RT*theta; forward AND reverse scale, preserving K."""
        if isinstance(step, bool) or not isinstance(step, int) or not 0 <= step < 18:
            raise ValueError("Step index must be an integer in [0,17]")
        theta = float(theta)
        if not math.isfinite(theta):
            raise ValueError("TS perturbation must be finite")
        try:
            factor = math.exp(theta)
        except OverflowError as error:
            raise ValueError("TS perturbation exceeds the numeric range") from error
        f, b = self.forward.copy(), self.reverse.copy()
        f[step] *= factor
        b[step] *= factor
        shifts = list(self.ts_perturbations_theta)
        shifts[step] += theta
        kind = self.evidence_kind
        if self.thermochemistry is not None:
            kind = "computed_TS_perturbation" if any(shifts) else "computed"
        return replace(self, forward=f, reverse=b, evidence_kind=kind,
                       ts_perturbations_theta=tuple(shifts))


class MasterKinetics:
    def __init__(self, rates: MasterRates, *, tbuoh_activity: float,
                 allow_software_fixture: bool = False):
        if rates.evidence_kind == "software_fixture" and not allow_software_fixture:
            raise EvidenceError("Synthetic software rates are not chemical predictions; opt in explicitly")
        if not math.isfinite(tbuoh_activity) or tbuoh_activity < 0:
            raise ValueError("tBuOH activity must be specified, finite and nonnegative")
        self.parameters = rates
        self.kf, self.kr = rates.forward, rates.reverse
        self.tbuoh_activity = float(tbuoh_activity)

    def _signature(self):
        digest = hashlib.sha256()
        for array in (self.kf, self.kr, np.array([self.parameters.temperature, self.tbuoh_activity])):
            digest.update(np.asarray(array, dtype="<f8").tobytes())
        digest.update(self.parameters.evidence_kind.encode("utf-8"))
        return digest.hexdigest()

    def rates(self, concentrations):
        C, CA, CHB, CH, CHI, CP, D, CO, X, A, B, N, H, I, P, W, Z, U = _finite_vector(concentrations, 18, "concentrations")
        f, b, a = self.kf, self.kr, self.tbuoh_activity
        return np.array([
            f[0]*C*A-b[0]*CA, f[1]*CA-b[1]*CHB,
            f[2]*CHB-b[2]*CH*B, f[3]*CH*I-b[3]*CHI,
            f[4]*CHI-b[4]*CP, f[5]*CP-b[5]*C*P,
            f[6]*B*N-b[6]*H, f[7]*H-b[7]*I*W,
            U*(f[8]*B*N-b[8]*H), U*(f[9]*H-b[9]*I*W),
            a*(f[10]*H-b[10]*I*W), a*a*(f[11]*H-b[11]*I*W),
            a*a*a*(f[12]*H-b[12]*I*W), f[13]*CH*CH-b[13]*D,
            f[14]*C*B, f[15]*CA*B, f[16]*C*U, f[17]*CA*U,
        ])

    def rate_jacobian(self, concentrations):
        """Hand-coded polynomial derivatives, including at zero concentration.

        Never uses r/C, log(C), automatic differentiation or finite differences.
        Rows are reaction channels; columns are the fixed 18 species.
        """
        C, CA, CHB, CH, CHI, CP, D, CO, X, A, B, N, H, I, P, W, Z, U = _finite_vector(concentrations, 18, "concentrations")
        f, b, a = self.kf, self.kr, self.tbuoh_activity
        d = np.zeros((18, 18))
        d[0, [0, 9, 1]] = [f[0]*A, f[0]*C, -b[0]]
        d[1, [1, 2]] = [f[1], -b[1]]
        d[2, [2, 3, 10]] = [f[2], -b[2]*B, -b[2]*CH]
        d[3, [3, 13, 4]] = [f[3]*I, f[3]*CH, -b[3]]
        d[4, [4, 5]] = [f[4], -b[4]]
        d[5, [5, 0, 14]] = [f[5], -b[5]*P, -b[5]*C]
        d[6, [10, 11, 12]] = [f[6]*N, f[6]*B, -b[6]]
        d[7, [12, 13, 15]] = [f[7], -b[7]*W, -b[7]*I]
        d[8, [10, 11, 12, 17]] = [U*f[8]*N, U*f[8]*B, -U*b[8], f[8]*B*N-b[8]*H]
        d[9, [12, 13, 15, 17]] = [U*f[9], -U*b[9]*W, -U*b[9]*I, f[9]*H-b[9]*I*W]
        for j, factor in ((10, a), (11, a*a), (12, a*a*a)):
            d[j, [12, 13, 15]] = [factor*f[j], -factor*b[j]*W, -factor*b[j]*I]
        d[13, [3, 6]] = [2*f[13]*CH, -b[13]]
        d[14, [0, 10]] = [f[14]*B, f[14]*C]
        d[15, [1, 10]] = [f[15]*B, f[15]*CA]
        d[16, [0, 17]] = [f[16]*U, f[16]*C]
        d[17, [1, 17]] = [f[17]*U, f[17]*CA]
        return d

    def rhs(self, t, concentrations):
        r = self.rates(concentrations)
        # The eighteen ODEs are explicit to make the material ledger reviewable.
        return np.array([
            -r[0]+r[5]-r[14]-r[16], r[0]-r[1]-r[15]-r[17],
            r[1]-r[2], r[2]-r[3]-2*r[13], r[3]-r[4], r[4]-r[5],
            r[13], r[14]+r[15], r[16]+r[17], -r[0]+r[15]+r[17],
            r[2]-r[6]-r[8]-r[14]-r[15], -r[6]-r[8],
            r[6]+r[8]-r[7]-r[9]-r[10]-r[11]-r[12],
            -r[3]+r[7]+r[9]+r[10]+r[11]+r[12],
            r[5], r[7]+r[9]+r[10]+r[11]+r[12], r[14]+r[15],
            -r[16]-r[17],
        ])

    def jacobian(self, t, concentrations):
        return S @ self.rate_jacobian(concentrations)

    def integrate(self, initial, *, end_time=3600.0, points=500,
                  rtol=1e-9, atol=1e-12):
        c0 = _finite_vector(initial, 18, "initial concentrations", nonnegative=True)
        if not math.isfinite(end_time) or end_time <= 0:
            raise ValueError("Positive finite end time required")
        if isinstance(points, bool) or not isinstance(points, int) or points < 2:
            raise ValueError("At least two integer output points required")
        if any(not math.isfinite(x) or x <= 0 for x in (rtol, atol)):
            raise ValueError("Integration tolerances must be positive and finite")
        metal_total = float(METAL @ c0)
        if metal_total <= 0:
            raise ValueError("Positive initial metal inventory required")
        sol = solve_ivp(self.rhs, (0., end_time), c0, method="Radau",
                        jac=self.jacobian, rtol=rtol, atol=atol, dense_output=True)
        if not sol.success:
            raise RuntimeError(f"Master ODE failed: {sol.message}")
        t = np.linspace(0., end_time, points)
        y = sol.sol(t)
        audit = np.concatenate([sol.y, y], axis=1)
        if not np.isfinite(audit).all():
            raise RuntimeError("Integration produced nonfinite concentrations")
        metal_error = float(np.max(np.abs(METAL @ audit-metal_total)))
        base_error = float(np.max(np.abs(BASE_EQUIVALENTS @ audit-BASE_EQUIVALENTS @ c0)))
        if metal_error >= 1e-10 or base_error >= 1e-10:
            raise RuntimeError(f"Inventory drift: metal={metal_error}, base={base_error} M")
        minimum = float(audit.min())
        if minimum < -1e-10:
            raise RuntimeError(f"Nonphysical negative concentration {minimum} M; no clipping applied")
        flux = np.array([self.rates(c) for c in y.T])
        if not np.isfinite(flux).all():
            raise RuntimeError("Integration produced nonfinite fluxes")
        elemental_error = None
        if self.parameters.thermochemistry is not None:
            labels, ledger = self.parameters.thermochemistry.elemental_charge_ledger()
            elemental_error = float(np.max(np.abs(ledger @ (audit-c0[:, None]))))
            # Atom inventories scale with ligand size; compare a per-row residual.
            residual = np.max(np.abs(ledger @ (audit-c0[:, None])), axis=1)
            if np.any(residual > 1e-10*np.maximum(1, np.max(np.abs(ledger), axis=1))):
                raise RuntimeError(f"Elemental/charge inventory drift: {dict(zip(labels, residual))}")
        return MasterTrajectory(t, y.T, flux, metal_total, metal_error, base_error,
                                minimum, self.parameters.evidence_kind, sol,
                                self._signature(), elemental_error)

    def control_coefficients(self, trajectory):
        """Integrate exact tangent equations on the accepted dense trajectory.

        theta_j = -delta G_TS,j / RT. Both rate directions scale by exp(theta).
        Z'=J Z+S diag(r), Z(0)=0. This is finite-time batch control, not a
        steady-state Campbell DRC claim. Undefined log derivatives remain NaN.
        """
        if trajectory.model_signature != self._signature():
            raise ValueError("Tangent analysis requires a trajectory from the same rates, temperature and reservoir activity")
        base = trajectory.solution

        def rhs(t, flat):
            c = base.sol(t)
            z = flat.reshape((18, 18), order="F")
            return (self.jacobian(t, c) @ z + S*self.rates(c)[None, :]).ravel(order="F")

        def jac(t, flat):
            return kron(eye(18, format="csc"), csr_matrix(self.jacobian(t, base.sol(t))), format="csc")

        sens = solve_ivp(rhs, (trajectory.time[0], trajectory.time[-1]),
                         np.zeros(324), method="Radau", jac=jac,
                         t_eval=trajectory.time, rtol=1e-9, atol=1e-12)
        if not sens.success:
            raise RuntimeError(f"Tangent ODE failed: {sens.message}")
        if not np.isfinite(sens.y).all():
            raise RuntimeError("Tangent ODE produced nonfinite sensitivities")
        z = np.array([column.reshape((18, 18), order="F") for column in sens.y.T])
        rate_control = np.full((len(trajectory.time), 18), np.nan)
        selectivity_control = rate_control.copy()
        yield_control = rate_control.copy()
        for q, (c, r, dz) in enumerate(zip(trajectory.concentrations, trajectory.fluxes, z)):
            dr = self.rate_jacobian(c) @ dz + np.diag(r)
            desired, byproduct = r[5], r[14]+r[15]
            if desired > 1e-20:
                rate_control[q] = dr[5]/desired
            if desired > 1e-20 and byproduct >= 0 and desired+byproduct > 1e-20:
                selectivity_control[q] = dr[5]/desired-(dr[5]+dr[14]+dr[15])/(desired+byproduct)
            product_formed = c[14]-trajectory.concentrations[0, 14]
            if product_formed > 1e-20:
                yield_control[q] = dz[14]/product_formed
        return {
            "kind": "finite_time_batch_TS_energy_control",
            "steady_state_campbell_drc_established": False,
            "theta_definition": "minus_delta_G_TS_over_RT; same factor in both directions",
            "reporting_flux_floor_M_per_s": 1e-20,
            "reporting_product_floor_M": 1e-20,
            "concentration_sensitivities": z,
            "rate_control": rate_control,
            "instantaneous_selectivity_control": selectivity_control,
            "cumulative_product_control": yield_control,
            "maximum_metal_sensitivity_residual_M": float(np.max(np.abs(np.einsum("i,tij->tj", METAL, z)))),
            "maximum_base_sensitivity_residual_M": float(np.max(np.abs(np.einsum("i,tij->tj", BASE_EQUIVALENTS, z)))),
        }


@dataclass
class MasterTrajectory:
    time: np.ndarray
    concentrations: np.ndarray
    fluxes: np.ndarray
    metal_total: float
    metal_error_M: float
    base_equivalent_error_M: float
    minimum_concentration_M: float
    evidence_kind: str
    solution: object
    model_signature: str = ""
    elemental_charge_error_M: float | None = None

    @property
    def tof_per_second(self):
        return self.fluxes[:, 5]/self.metal_total

    @property
    def active_metal_fraction(self):
        return self.concentrations[:, :6].sum(axis=1)/self.metal_total

    @property
    def cleavage_fragment_formed_M(self):
        return self.concentrations[:, 8]-self.concentrations[0, 8]


def software_fixture_rates(temperature=383.15):
    """Deterministic stiffness/conservation fixture; NOT calculated chemistry.

    Arbitrary activation enthalpies only exercise temperature plumbing. No
    experimental fit, mechanistic ranking or catalyst identity is implied.
    """
    temperature = float(temperature)
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("Fixture temperature must be positive and finite")
    f = np.array([30, 2, 10, 25, 1.5, 8, 1, .08, 3, 2, .3, .4, .2, 40, .004, .002, .003, .002])
    b = np.array([3, .1, 2, 2, .1, .5, .08, .01, .24, .25, .0375, .05, .025, .002, 0, 0, 0, 0])
    # Reverse ratios of parallel condensation channels match deliberately.
    # This is a solver fixture, not a thermochemical certificate.
    activation = np.linspace(5., 12., 18)
    factors = np.exp(-activation/R_KCAL*(1/float(temperature)-1/383.15))
    return MasterRates(f*factors, b*factors, float(temperature), "software_fixture")


def software_fixture_initial(*, alcohol_M=1.0, catalyst_M=.01, free_base_equiv=.05):
    values = (alcohol_M, catalyst_M, free_base_equiv)
    if any(not math.isfinite(x) or x < 0 for x in values) or alcohol_M <= 0 or catalyst_M <= 0:
        raise ValueError("Fixture concentrations/equivalents must be finite and nonnegative, with alcohol/catalyst positive")
    c = np.zeros(18)
    c[0], c[9], c[11], c[17] = catalyst_M, alcohol_M, alcohol_M, alcohol_M*free_base_equiv
    return c
