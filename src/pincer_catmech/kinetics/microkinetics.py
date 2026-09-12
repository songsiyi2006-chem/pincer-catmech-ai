"""Evidence-gated borrowing-hydrogen mass action in mol/L and seconds.

Six reversible steps contain eleven species. Catalyst, benzyl-fragment and
aniline-fragment balances eliminate catalyst, alcohol and amine: eight ODEs.
Ligand-bound nitrogen/carbon are constant with catalyst inventory; the two
organic-fragment balances therefore also conserve total organic C and N.
All free energies are computed 1 M qRRHO values at the simulated temperature.
No missing rate, barrier, temperature dependence, or base order is inferred.

Batch barrier sensitivities are explicitly finite-time apparent DRC, not a
claim of steady-state Campbell DRC (doi:10.1021/acscatal.7b00115).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from collections import Counter
import hashlib
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.integrate import solve_ivp

R_KCAL = 8.31446261815324 / 4184.0
KB_OVER_H = 1.380649e-23 / 6.62607015e-34
SPECIES = ("alcohol", "cat", "complex_1", "aldehyde", "cat_h2", "amine",
           "hemiaminal", "imine", "water", "complex_2", "product")
INDEX = {name: i for i, name in enumerate(SPECIES)}
DYNAMIC = tuple(i for i in range(len(SPECIES)) if i not in (0, 1, 5))
STEPS = {
    "association_1": ({"alcohol": 1, "cat": 1}, {"complex_1": 1}),
    "dehydrogenation": ({"complex_1": 1}, {"aldehyde": 1, "cat_h2": 1}),
    "hemiaminal_formation": ({"aldehyde": 1, "amine": 1}, {"hemiaminal": 1}),
    "dehydration": ({"hemiaminal": 1}, {"imine": 1, "water": 1}),
    "association_2": ({"imine": 1, "cat_h2": 1}, {"complex_2": 1}),
    "hydrogenation": ({"complex_2": 1}, {"product": 1, "cat": 1}),
}
BALANCES = np.array([
    [0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0],
    [1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1],
    [0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1],
], dtype=float)


class EvidenceError(ValueError):
    """Required accepted computed thermochemistry is absent or inconsistent."""


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


@dataclass(frozen=True)
class ComputedFreeEnergy:
    """Traceable scalar. Flags attest supplied calculation review, not new QM."""
    gibbs_kcal_mol: float
    temperature: float
    source_path: str
    source_sha256: str
    method: str
    accepted: bool
    stationarity_verified: bool
    imaginary_frequencies_cm1: tuple[float, ...] = ()
    connectivity_verified: bool = False
    evidence_kind: str = "computed"
    standard_state: str = "1M"
    solvent: str | None = None
    solvation_state: str | None = None
    chemical_identity_verified: bool = False
    reaction_mode_verified: bool = False
    composition: tuple[tuple[str, int], ...] = ()
    charge: int | None = None

    @property
    def protocol(self) -> tuple[str, str | None, str | None]:
        return self.method, self.solvent, self.solvation_state

    def validate(self, temperature: float, *, transition_state: bool = False,
                 dehydrogenation: bool = False) -> None:
        if not math.isfinite(temperature) or temperature <= 0 or not math.isfinite(self.temperature) or self.temperature <= 0:
            raise EvidenceError("Computed free-energy temperatures must be positive and finite")
        if self.evidence_kind != "computed" or self.accepted is not True or self.stationarity_verified is not True:
            raise EvidenceError("Only accepted stationary computed free energies are permitted")
        if not self.method.strip() or self.standard_state != "1M":
            raise EvidenceError("Method and 1 M free-energy reference are required")
        if not math.isfinite(self.gibbs_kcal_mol) or not math.isclose(
            self.temperature, temperature, abs_tol=1e-6, rel_tol=0
        ):
            raise EvidenceError("Free energy must be finite and computed at the requested temperature")
        source = Path(self.source_path)
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != self.source_sha256:
            raise EvidenceError("Computed source artifact is absent or its SHA256 differs")
        if not self.solvent or not self.solvation_state:
            raise EvidenceError("Explicit solvent and solvation reference state are required; use gas/none for a gas model")
        if self.chemical_identity_verified is not True:
            raise EvidenceError("The full intended chemical identity must be independently verified")
        if not self.composition or len(dict(self.composition)) != len(self.composition) or any(
            not isinstance(symbol, str) or not symbol.isalpha() or
            isinstance(count, bool) or not isinstance(count, int) or count < 1
            for symbol, count in self.composition
        ) or isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise EvidenceError("Explicit elemental composition and integral molecular charge are required")
        frequencies = self.imaginary_frequencies_cm1
        if any(not math.isfinite(v) or v >= 0 for v in frequencies):
            raise EvidenceError("Imaginary-frequency metadata must contain finite negative modes")
        if transition_state:
            if len(frequencies) != 1 or self.connectivity_verified is not True or self.reaction_mode_verified is not True:
                raise EvidenceError("Transition state needs exactly one imaginary mode, verified reaction motion and connectivity")
            if dehydrogenation and not -1800 <= frequencies[0] <= -300:
                raise EvidenceError("Dehydrogenation mode fails the required [-1800,-300] cm^-1 criterion")
        elif frequencies:
            raise EvidenceError("Intermediate has imaginary modes")


@dataclass(frozen=True)
class ThermodynamicSnapshot:
    temperature: float
    species: Mapping[str, ComputedFreeEnergy]
    transition_states: Mapping[str, ComputedFreeEnergy]
    base_assisted_transition_states: Mapping[str, ComputedFreeEnergy] | None = None
    base_free_energy: ComputedFreeEnergy | None = None
    co_shuttle_free_energy: ComputedFreeEnergy | None = None

    def validate(self, require_base_path: bool = False) -> None:
        temperature = _positive(self.temperature, "temperature")
        if set(self.species) != set(SPECIES):
            raise EvidenceError(f"All eleven species are required; missing {sorted(set(SPECIES)-set(self.species))}")
        if set(self.transition_states) != set(STEPS):
            raise EvidenceError(f"All six step barriers are required; missing {sorted(set(STEPS)-set(self.transition_states))}")
        for energy in self.species.values():
            energy.validate(temperature)
        for step, energy in self.transition_states.items():
            energy.validate(temperature, transition_state=True, dehydrogenation=step == "dehydrogenation")
        paths = self.base_assisted_transition_states or {}
        records = list(self.species.values()) + list(self.transition_states.values()) + list(paths.values())
        if self.base_free_energy is not None:
            records.append(self.base_free_energy)
        if self.co_shuttle_free_energy is not None:
            records.append(self.co_shuttle_free_energy)
        if len({record.protocol for record in records}) != 1:
            raise EvidenceError("All network energies must use one computational method, solvent and reference-state protocol")
        if require_base_path and not paths:
            raise EvidenceError("Base grid requires a computed base-assisted path; no kinetic order is assumed")
        if not paths and (self.base_free_energy is not None or self.co_shuttle_free_energy is not None):
            raise EvidenceError("Shuttle reference energies require a declared base-assisted path")
        if set(paths) - {"hemiaminal_formation", "dehydration"}:
            raise EvidenceError("Only explicitly modeled single-base condensation paths are supported")
        if paths:
            if self.base_free_energy is None:
                raise EvidenceError("A computed base reference energy is required")
            self.base_free_energy.validate(temperature)
            if self.co_shuttle_free_energy is not None:
                self.co_shuttle_free_energy.validate(temperature)
            for energy in paths.values():
                energy.validate(temperature, transition_state=True)
        for step, (reactants, products) in STEPS.items():
            references = [[(self.species[name], n) for name, n in side.items()]
                          for side in (reactants, products)]
            expected = _composition_balance(references[0])
            if expected != _composition_balance(references[1]):
                raise EvidenceError(f"{step} does not conserve elemental composition and charge")
            if expected != _composition_balance([(self.transition_states[step], 1)]):
                raise EvidenceError(f"{step} transition-state composition differs from its reactants")
            if step in paths:
                shuttle = [(self.base_free_energy, 1)]
                if self.co_shuttle_free_energy is not None:
                    shuttle.append((self.co_shuttle_free_energy, 1))
                if _composition_balance(references[0] + shuttle) != _composition_balance([(paths[step], 1)]):
                    raise EvidenceError(f"{step} base-assisted TS does not contain the declared shuttle composition/charge")


def _composition_balance(records: Sequence[tuple[ComputedFreeEnergy, int]]) -> tuple[dict[str, int], int]:
    elements: Counter[str] = Counter()
    charge = 0
    for record, coefficient in records:
        elements.update({symbol: count*coefficient for symbol, count in record.composition})
        charge += record.charge*coefficient
    return dict(elements), charge


def eyring_rate_constant(barrier_kcal_mol: float, temperature: float, molecularity: int,
                         standard_concentration_mol_l: float = 1.0) -> float:
    """k has units M**(1-molecularity)/s; transmission coefficient is one."""
    temperature = _positive(temperature, "temperature")
    c0 = _positive(standard_concentration_mol_l, "standard concentration")
    if isinstance(molecularity, bool) or not isinstance(molecularity, int) or molecularity < 1:
        raise ValueError("molecularity must be a positive integer")
    if not math.isfinite(barrier_kcal_mol) or barrier_kcal_mol < 0:
        raise ValueError("Barrier must be finite and nonnegative; barrierless capture needs another model")
    log_k = math.log(KB_OVER_H * temperature) - barrier_kcal_mol / (R_KCAL * temperature)
    log_k += (1 - molecularity) * math.log(c0)
    if log_k > math.log(np.finfo(float).max):
        raise ValueError("Rate constant exceeds floating-point range")
    return math.exp(log_k)


@dataclass(frozen=True)
class ReactionChannel:
    name: str
    reactants: tuple[tuple[int, int], ...]
    products: tuple[tuple[int, int], ...]
    k_forward: float
    k_reverse: float
    base_activity: float = 1.0

    def __post_init__(self) -> None:
        if not self.name or not self.reactants or not self.products:
            raise ValueError("Named channel and both stoichiometric sides are required")
        for side in (self.reactants, self.products):
            if len({i for i, _ in side}) != len(side):
                raise ValueError("Duplicate species in stoichiometric side")
            if any(not isinstance(i, int) or not 0 <= i < 11 or
                   not isinstance(n, int) or n <= 0 for i, n in side):
                raise ValueError("Invalid stoichiometric index or coefficient")
        if any(not math.isfinite(v) or v < 0 for v in (self.k_forward, self.k_reverse, self.base_activity)):
            raise ValueError("Rates and base activity must be finite and nonnegative")

    def flux(self, concentrations: np.ndarray) -> float:
        forward = self.k_forward * math.prod(max(0.0, concentrations[i])**n for i, n in self.reactants)
        reverse = self.k_reverse * math.prod(max(0.0, concentrations[i])**n for i, n in self.products)
        return self.base_activity * (forward - reverse)


@dataclass(frozen=True)
class MicrokineticModel:
    channels: tuple[ReactionChannel, ...]
    temperature: float
    evidence_level: str
    source_snapshot: ThermodynamicSnapshot | None = None
    free_base_concentration_mol_l: float | None = None
    co_shuttle_concentration_mol_l: float | None = None

    def __post_init__(self) -> None:
        _positive(self.temperature, "temperature")
        if not self.channels or len({c.name for c in self.channels}) != len(self.channels):
            raise ValueError("At least one uniquely named channel is required")
        if self.evidence_level not in ("computed_accepted_thermochemistry", "analytical_unit_test"):
            raise ValueError("Explicit computed or analytical-test evidence level required")
        if self.evidence_level == "computed_accepted_thermochemistry":
            if self.source_snapshot is None:
                raise EvidenceError("Computed rates require their complete accepted source snapshot")
            self.source_snapshot.validate()
            if not math.isclose(self.source_snapshot.temperature, self.temperature, rel_tol=0, abs_tol=1e-6):
                raise EvidenceError("Model and source snapshot temperatures differ")
            expected_names = set(STEPS) | {step+":base" for step in self.source_snapshot.base_assisted_transition_states or {}}
            if {channel.name for channel in self.channels} != expected_names:
                raise EvidenceError("Computed model channels must contain every source-snapshot step")
            for channel in self.channels:
                reactants, products = STEPS[channel.name.split(":", 1)[0]]
                if dict(channel.reactants) != {INDEX[name]: n for name, n in reactants.items()} or dict(channel.products) != {
                    INDEX[name]: n for name, n in products.items()
                }:
                    raise EvidenceError("Computed channel stoichiometry differs from its composition-verified source step")
        if not np.allclose(BALANCES @ self.stoichiometry, 0, atol=0, rtol=0):
            raise ValueError("Network violates catalyst/carbon-fragment/nitrogen-fragment conservation")

    @property
    def stoichiometry(self) -> np.ndarray:
        matrix = np.zeros((11, len(self.channels)))
        for j, channel in enumerate(self.channels):
            for i, n in channel.reactants:
                matrix[i, j] -= n
            for i, n in channel.products:
                matrix[i, j] += n
        return matrix


def build_model(snapshot: ThermodynamicSnapshot, *, base_concentration_mol_l: float = 0,
                require_base_path: bool = False,
                co_shuttle_concentration_mol_l: float | None = None) -> MicrokineticModel:
    """Build reversible rates from shared TS and intermediate absolute free energies."""
    snapshot.validate(require_base_path=require_base_path)
    if not math.isfinite(base_concentration_mol_l) or base_concentration_mol_l < 0:
        raise ValueError("Base concentration must be finite and nonnegative")
    if base_concentration_mol_l > 0 and not snapshot.base_assisted_transition_states:
        raise EvidenceError("A nonzero free-base input requires an explicitly computed base-assisted path")
    if co_shuttle_concentration_mol_l is not None and snapshot.co_shuttle_free_energy is None:
        raise EvidenceError("A co-shuttle concentration requires its computed reference and explicit path model")
    if snapshot.co_shuttle_free_energy is not None:
        if co_shuttle_concentration_mol_l is None:
            raise EvidenceError("A two-species proton shuttle requires its separately specified co-shuttle concentration")
        if not math.isfinite(co_shuttle_concentration_mol_l) or co_shuttle_concentration_mol_l < 0:
            raise ValueError("Co-shuttle concentration must be finite and nonnegative")
    channels = []
    for step, (reactants, products) in STEPS.items():
        paths = [(step, snapshot.transition_states[step], False)]
        if step in (snapshot.base_assisted_transition_states or {}):
            paths.append((step + ":base", snapshot.base_assisted_transition_states[step], True))
        for name, transition_state, assisted in paths:
            base_g = snapshot.base_free_energy.gibbs_kcal_mol if assisted else 0.0
            activity = base_concentration_mol_l / 1.0 if assisted else 1.0
            if assisted and snapshot.co_shuttle_free_energy is not None:
                base_g += snapshot.co_shuttle_free_energy.gibbs_kcal_mol
                activity *= co_shuttle_concentration_mol_l / 1.0
            barriers = [transition_state.gibbs_kcal_mol - base_g - sum(
                n * snapshot.species[species].gibbs_kcal_mol for species, n in side.items()
            ) for side in (reactants, products)]
            rates = [eyring_rate_constant(barrier, snapshot.temperature, sum(side.values()))
                     for barrier, side in zip(barriers, (reactants, products))]
            channels.append(ReactionChannel(
                name, tuple((INDEX[s], n) for s, n in reactants.items()),
                tuple((INDEX[s], n) for s, n in products.items()), rates[0], rates[1],
                activity,
            ))
    return MicrokineticModel(tuple(channels), snapshot.temperature, "computed_accepted_thermochemistry",
                            snapshot, base_concentration_mol_l, co_shuttle_concentration_mol_l)


def initial_concentrations(*, alcohol_mol_l: float = 1.0, amine_mol_l: float = 1.0,
                           catalyst_loading_mol_percent: float = 1.0) -> dict[str, float]:
    alcohol = _positive(alcohol_mol_l, "alcohol concentration")
    amine = _positive(amine_mol_l, "amine concentration")
    loading = _positive(catalyst_loading_mol_percent, "catalyst loading")
    state = dict.fromkeys(SPECIES, 0.0)
    state.update(alcohol=alcohol, amine=amine, cat=alcohol * loading / 100.0)
    return state


def _expand(reduced: np.ndarray, balances: np.ndarray) -> np.ndarray:
    full = np.zeros(11)
    full[list(DYNAMIC)] = reduced
    full[1] = balances[0] - full[2] - full[4] - full[9]
    full[0] = balances[1] - sum(full[i] for i in (2, 3, 6, 7, 9, 10))
    full[5] = balances[2] - sum(full[i] for i in (6, 7, 9, 10))
    return full


@dataclass(frozen=True)
class KineticResult:
    time_seconds: np.ndarray
    concentrations_mol_l: np.ndarray
    average_tof_per_second: float
    instantaneous_tof_per_second: float
    conservation_max_abs_mol_l: float
    solver_method: str
    rhs_evaluations: int
    evidence_level: str
    integrated_species: tuple[str, ...] = tuple(SPECIES[i] for i in DYNAMIC)
    temperature: float | None = None
    initial_concentrations_mol_l: tuple[float, ...] = ()
    free_base_concentration_mol_l: float | None = None
    co_shuttle_concentration_mol_l: float | None = None
    computational_protocol: tuple[str, str | None, str | None] | None = None


def solve_batch(model: MicrokineticModel, initial: Mapping[str, float], *,
                end_time_seconds: float, n_output: int = 101, method: str = "Radau",
                rtol: float = 1e-8, atol: float = 1e-11) -> KineticResult:
    """Closed batch eight-ODE integration; output reconstructs all eleven species."""
    if set(initial) != set(SPECIES):
        raise ValueError("Initial state must specify all eleven species")
    y0 = np.array([initial[s] for s in SPECIES], dtype=float)
    if np.any(~np.isfinite(y0)) or np.any(y0 < 0):
        raise ValueError("Initial concentrations must be finite and nonnegative")
    if method not in ("Radau", "BDF") or isinstance(n_output, bool) or not isinstance(n_output, int) or n_output < 2:
        raise ValueError("Use Radau/BDF and at least two output points")
    end_time = _positive(end_time_seconds, "end time")
    _positive(rtol, "rtol")
    _positive(atol, "atol")
    inventory = BALANCES @ y0
    _positive(inventory[0], "total catalyst")
    matrix = model.stoichiometry

    def rhs(_time: float, reduced: np.ndarray) -> np.ndarray:
        concentrations = _expand(reduced, inventory)
        return (matrix @ np.array([c.flux(concentrations) for c in model.channels]))[list(DYNAMIC)]

    times = np.linspace(0, end_time, n_output)
    solution = solve_ivp(rhs, (0, end_time), y0[list(DYNAMIC)], method=method,
                         t_eval=times, rtol=rtol, atol=atol)
    if not solution.success or solution.y.shape[1] != n_output:
        raise RuntimeError(f"Stiff integration failed: {solution.message}")
    concentrations = np.stack([_expand(column, inventory) for column in solution.y.T])
    if np.any(~np.isfinite(concentrations)) or concentrations.min() < -max(100*atol, 10*rtol*y0.max()):
        raise RuntimeError("Integration produced materially negative or nonfinite concentrations")
    final_derivative = matrix @ np.array([c.flux(concentrations[-1]) for c in model.channels])
    deviation = float(np.max(np.abs(concentrations @ BALANCES.T - inventory)))
    return KineticResult(times, concentrations,
                         float((concentrations[-1, 10] - y0[10]) / (end_time * inventory[0])),
                         float(final_derivative[10] / inventory[0]), deviation,
                         method, solution.nfev, model.evidence_level,
                         temperature=model.temperature, initial_concentrations_mol_l=tuple(y0),
                         free_base_concentration_mol_l=model.free_base_concentration_mol_l,
                         co_shuttle_concentration_mol_l=model.co_shuttle_concentration_mol_l,
                         computational_protocol=(next(iter(model.source_snapshot.species.values())).protocol
                                                 if model.source_snapshot is not None else None))


@dataclass(frozen=True)
class RateSensitivity:
    channel: str
    value: float | None
    status: str
    perturbation_kcal_mol: float
    definition: str = "finite_time_apparent_DRC_of_average_TOF"


def campbell_drc(model: MicrokineticModel, initial: Mapping[str, float], *,
                  end_time_seconds: float, perturbation_kcal_mol: float = 0.02,
                  method: str = "Radau") -> tuple[RateSensitivity, ...]:
    """-RT*d ln(average TOF)/dG_TS; both directions shifted, intermediates fixed.

    Finite batch time gives apparent DRC, not steady-state Campbell DRC.
    No sum-to-one identity is imposed on this finite-time sensitivity.
    """
    delta = _positive(perturbation_kcal_mol, "barrier perturbation")
    results = []
    for index, channel in enumerate(model.channels):
        rates = []
        for sign in (-1, 1):
            factor = math.exp(-sign * delta / (R_KCAL * model.temperature))
            changed = list(model.channels)
            changed[index] = replace(channel, k_forward=channel.k_forward*factor,
                                     k_reverse=channel.k_reverse*factor)
            result = solve_batch(replace(model, channels=tuple(changed)), initial,
                                 end_time_seconds=end_time_seconds, n_output=2, method=method)
            rates.append(result.average_tof_per_second)
        if min(rates) <= 0:
            results.append(RateSensitivity(channel.name, None, "undefined_nonpositive_net_TOF", delta))
        else:
            sensitivity = -R_KCAL*model.temperature*(math.log(rates[1])-math.log(rates[0]))/(2*delta)
            results.append(RateSensitivity(channel.name, sensitivity, "computed", delta))
    return tuple(results)


def campaign_readiness(snapshots: Mapping[float, ThermodynamicSnapshot], *,
                       temperatures: Sequence[float] = tuple(np.linspace(340, 440, 11)),
                       require_base_path: bool = True) -> dict:
    if not len(temperatures):
        raise ValueError("At least one temperature is required")
    problems = []
    for temperature in temperatures:
        temperature = _positive(temperature, "temperature")
        if temperature not in snapshots:
            problems.append(f"Missing complete accepted thermochemistry at {temperature:g} K")
        else:
            try:
                if not math.isclose(snapshots[temperature].temperature, temperature, abs_tol=1e-6, rel_tol=0):
                    raise EvidenceError("Snapshot key and temperature differ")
                snapshots[temperature].validate(require_base_path)
            except (EvidenceError, ValueError) as error:
                problems.append(f"{temperature:g} K: {error}")
    return {"status": "ready" if not problems else "insufficient_computed_data", "reasons": problems}


def run_condition_grid(snapshots: Mapping[float, ThermodynamicSnapshot], *, end_time_seconds: float,
                        temperatures: Sequence[float] = tuple(np.linspace(340, 440, 11)),
                        loadings_mol_percent: Sequence[float] = tuple(np.linspace(0.1, 2.0, 5)),
                        base_equivalents: Sequence[float] = (0.01, 0.05, 0.10, 0.20),
                        alcohol_mol_l: float = 1.0, amine_mol_l: float = 1.0,
                        calculate_drc: bool = True, method: str = "Radau",
                        free_base_fraction: float | None = None,
                        co_shuttle_concentration_mol_l: float | None = None) -> list[dict]:
    """Explicit study conditions; no execution when any required evidence is absent."""
    readiness = campaign_readiness(snapshots, temperatures=temperatures, require_base_path=True)
    if readiness["status"] != "ready":
        raise EvidenceError("; ".join(readiness["reasons"]))
    if free_base_fraction is None:
        raise EvidenceError("Total base equivalents require an explicit free-base fraction/speciation assumption")
    if not math.isfinite(free_base_fraction) or not 0 <= free_base_fraction <= 1:
        raise ValueError("Free-base fraction must be between zero and one")
    if not len(temperatures) or not len(loadings_mol_percent) or not len(base_equivalents):
        raise ValueError("Condition axes must be nonempty")
    rows = []
    for temperature in temperatures:
        for loading in loadings_mol_percent:
            initial = initial_concentrations(alcohol_mol_l=alcohol_mol_l, amine_mol_l=amine_mol_l,
                                              catalyst_loading_mol_percent=loading)
            for base in base_equivalents:
                base = _positive(base, "base equivalents")
                model = build_model(snapshots[temperature], base_concentration_mol_l=base*alcohol_mol_l*free_base_fraction,
                                    require_base_path=True, co_shuttle_concentration_mol_l=co_shuttle_concentration_mol_l)
                result = solve_batch(model, initial, end_time_seconds=end_time_seconds, method=method)
                sensitivities = campbell_drc(model, initial, end_time_seconds=end_time_seconds, method=method) if calculate_drc else ()
                rows.append({"temperature": temperature, "loading_mol_percent": loading,
                             "base_equivalents": base, "free_base_fraction": free_base_fraction,
                             "total_base_concentration_mol_l": base*alcohol_mol_l,
                             "free_base_concentration_mol_l": base*alcohol_mol_l*free_base_fraction,
                             "co_shuttle_concentration_mol_l": co_shuttle_concentration_mol_l,
                             "base_model": "fixed free-species activities; no ion-pair or activation equilibrium inferred",
                             "result": result, "drc": sensitivities})
    return rows


__all__ = ["SPECIES", "STEPS", "BALANCES", "ComputedFreeEnergy", "ThermodynamicSnapshot",
           "EvidenceError", "eyring_rate_constant", "ReactionChannel", "MicrokineticModel",
           "build_model", "initial_concentrations", "solve_batch", "campbell_drc",
           "campaign_readiness", "run_condition_grid", "KineticResult", "RateSensitivity"]
