"""Thermodynamically consistent reversible mass action; no physical rate certificate.

Energies: kcal/mol. Concentrations: mol/L. Activities: dimensionless relative
to a common c0 in mol/L. Fluxes: mol/L/s. Only ideal closed-system free energy
has a Lyapunov diagnostic here; arbitrary activity coefficients do not.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Mapping
import math

import numpy as np

R_KCAL = 8.31446261815324 / 4184.0
K_B = 1.380649e-23
H_PLANCK = 6.62607015e-34
ELEMENT_SYMBOLS = set(('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
                      'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La '
                      'Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi '
                      'Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt '
                      'Ds Rg Cn Nh Fl Mc Lv Ts Og').split())


def _scalar(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be a finite real scalar')
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real scalar') from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f'{name} must be finite' + (' and positive' if positive else ''))
    return result


def _array(value, shape, name):
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real array') from exc
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f'{name} must have shape {shape} and finite entries')
    # Immutable bytes back the view; callers cannot re-enable WRITEABLE.
    return np.frombuffer(array.tobytes(), dtype=float).reshape(shape)


def _integers(value, shape, name, *, nonnegative=False):
    raw = np.asarray(value, dtype=object)
    if raw.shape != shape or any(isinstance(x, (bool, np.bool_)) or
            not isinstance(x, (int, np.integer)) for x in raw.flat):
        raise ValueError(f'{name} must have shape {shape} and integer entries')
    if any(abs(int(x)) > 10**6 or (nonnegative and int(x) < 0) for x in raw.flat):
        raise ValueError(f'{name} has negative or unsupported integer entries')
    array = np.asarray(raw, dtype=np.int64)
    return np.frombuffer(array.tobytes(), dtype=np.int64).reshape(shape)


def _names(values, name):
    result = tuple(values)
    if not result or any(not isinstance(x, str) or not x.strip() for x in result) or len(set(result)) != len(result):
        raise ValueError(f'{name} must be nonempty unique names')
    return result


def _finite_result(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f'{name} exceeds finite floating-point range')
    return value


def _exp(log_values):
    try:
        with np.errstate(over='raise', invalid='raise', under='ignore'):
            result = np.exp(log_values)
    except FloatingPointError as exc:
        raise ValueError('Rate exceeds finite floating-point range') from exc
    return _finite_result(result, 'Rate')


@dataclass(frozen=True)
class ThermodynamicNetwork:
    """General mathematical network, not an experimentally calibrated catalyst.

    alpha and beta have shape (species, reactions). A single common transition
    standard chemical potential is used for both directions of each reaction.
    Compositions are explicit elemental inventories; every listed species must
    contain at least one element. No implicit proton, counterion, or chemostat
    is silently introduced. Shared energies establish thermodynamic consistency,
    not the validity of transition-state theory for a particular chemistry.
    """

    species: tuple[str, ...]
    reactions: tuple[str, ...]
    species_mu0: np.ndarray
    transition_mu0: np.ndarray
    alpha: np.ndarray
    beta: np.ndarray
    compositions: tuple[Mapping[str, int], ...]
    charges: np.ndarray
    temperature_K: float
    standard_concentration_M: float = 1.0
    physical_rates_validated: bool = field(default=False, init=False)
    _elements: tuple[str, ...] = field(init=False, repr=False)
    _inventory: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        species, reactions = _names(self.species, 'species'), _names(self.reactions, 'reactions')
        ns, nr = len(species), len(reactions)
        alpha = _integers(self.alpha, (ns, nr), 'alpha', nonnegative=True)
        beta = _integers(self.beta, (ns, nr), 'beta', nonnegative=True)
        charges = _integers(self.charges, (ns,), 'charges')
        if len(self.compositions) != ns:
            raise ValueError('Every species requires an explicit elemental composition')
        compositions = []
        for comp in self.compositions:
            if not isinstance(comp, Mapping) or not comp:
                raise ValueError('Species elemental inventory cannot be empty')
            copied = {}
            for element, count in comp.items():
                if not isinstance(element, str) or element not in ELEMENT_SYMBOLS:
                    raise ValueError('Element labels must be explicit case-sensitive symbols')
                val = _integers([count], (1,), 'element count', nonnegative=True)[0]
                if val <= 0:
                    raise ValueError('Listed element counts must be positive')
                copied[element] = int(val)
            compositions.append(MappingProxyType(copied))
        elements = tuple(sorted({element for comp in compositions for element in comp}))
        inventory = np.asarray([[comp.get(e, 0) for comp in compositions] for e in elements], dtype=np.int64)
        nu = beta - alpha
        # Object arithmetic avoids integer overflow in conservation assertions.
        if np.any(inventory.astype(object) @ nu.astype(object)):
            raise ValueError('Reaction violates elemental conservation')
        if np.any(charges.astype(object) @ nu.astype(object)):
            raise ValueError('Reaction violates charge conservation')
        if np.any(np.all(nu == 0, axis=0)):
            raise ValueError('A reaction must change at least one species')
        temperature = _scalar(self.temperature_K, 'temperature_K', positive=True)
        concentration = _scalar(self.standard_concentration_M, 'standard_concentration_M', positive=True)
        for name, value in dict(species=species, reactions=reactions, alpha=alpha, beta=beta,
                species_mu0=_array(self.species_mu0, (ns,), 'species_mu0'),
                transition_mu0=_array(self.transition_mu0, (nr,), 'transition_mu0'), charges=charges,
                compositions=tuple(compositions), temperature_K=temperature,
                standard_concentration_M=concentration, _elements=elements,
                _inventory=_integers(inventory, inventory.shape, 'inventory', nonnegative=True)).items():
            object.__setattr__(self, name, value)

    @property
    def nu(self):
        return self.beta - self.alpha

    @property
    def rt(self):
        return R_KCAL * self.temperature_K

    def _activities(self, activities, *, interior=False):
        a = _array(activities, (len(self.species),), 'activities')
        if np.any(a < 0) or (interior and np.any(a <= 0)):
            raise ValueError('Activities must be ' + ('strictly positive in the interior' if interior else 'nonnegative'))
        return a

    def _concentrations(self, concentrations, *, interior=False):
        c = _array(concentrations, (len(self.species),), 'concentrations_M')
        if np.any(c < 0) or (interior and np.any(c <= 0)):
            raise ValueError('Concentrations must be ' + ('strictly positive in the interior' if interior else 'nonnegative'))
        return c

    def log_frequency_constants(self):
        """Forward/reverse log constants, each s^-1 for dimensionless activities."""
        prefactor = math.log(K_B / H_PLANCK) + math.log(self.temperature_K)
        with np.errstate(over='ignore', invalid='ignore'):
            forward = prefactor - (self.transition_mu0 - self.alpha.T @ self.species_mu0) / self.rt
            reverse = prefactor - (self.transition_mu0 - self.beta.T @ self.species_mu0) / self.rt
        return (_finite_result(forward, 'forward log constants'),
                _finite_result(reverse, 'reverse log constants'))

    def log_equilibrium_constants(self):
        """Dimensionless K for prod(a_i ** nu_i); shared species standard state."""
        with np.errstate(over='ignore', invalid='ignore'):
            return _finite_result(-(self.nu.T @ self.species_mu0) / self.rt, 'log K')

    def concentration_rate_constants(self):
        """Conventional k multiplying products of molar concentrations.

        Units differ by molecularity: M^(1-sum(alpha))*s^-1 and the corresponding
        beta order. They must not be treated as uniformly first-order constants.
        """
        f, r = self.log_frequency_constants()
        pf, pr = 1 - self.alpha.sum(axis=0), 1 - self.beta.sum(axis=0)
        return dict(forward=_exp(f + pf * math.log(self.standard_concentration_M)),
                    reverse=_exp(r + pr * math.log(self.standard_concentration_M)),
                    forward_molarity_power=pf, reverse_molarity_power=pr, time_power=-1)

    def log_unidirectional_fluxes(self, activities):
        """Log(M/s); a zero reactant activity gives -inf, with 0**0 = 1."""
        a = self._activities(activities)
        # A common transition-state offset reduces cancellation near equilibrium.
        prefactor = math.log(K_B / H_PLANCK) + math.log(self.temperature_K)
        common = _finite_result(prefactor + math.log(self.standard_concentration_M)
                                - self.transition_mu0 / self.rt, 'common flux offset')
        logs = []
        for stoich in (self.alpha, self.beta):
            value = common.copy()
            for j in range(len(self.reactions)):
                used = stoich[:, j] > 0
                if np.any(a[used] == 0):
                    value[j] = -math.inf
                else:
                    value[j] += np.dot(stoich[used, j], self.species_mu0[used] / self.rt + np.log(a[used]))
            if np.any(np.isnan(value)) or np.any(np.isposinf(value)):
                raise ValueError('Flux logarithm exceeds finite floating-point range')
            logs.append(value)
        return tuple(logs)

    def unidirectional_fluxes(self, activities):
        return tuple(_exp(x) for x in self.log_unidirectional_fluxes(activities))

    def flux_from_activities(self, activities):
        """Net reversible flux M/s. Supplied nonideal activities confer no Lyapunov proof."""
        lf, lr = self.log_unidirectional_fluxes(activities)
        result = np.zeros(len(self.reactions))
        for j, (f, r) in enumerate(zip(lf, lr)):
            if f == r:  # Includes both -inf at a zero-activity boundary.
                continue
            if f > r:
                result[j] = float(_exp(f)) * (-math.expm1(r - f))
            else:
                result[j] = -float(_exp(r)) * (-math.expm1(f - r))
        return _finite_result(result, 'net flux')

    def rhs(self, concentrations_M):
        """Ideal closed-network vector field M/s, including nonnegative boundary.

        This function is not an ODE integrator and does not clip negative states.
        """
        c = self._concentrations(concentrations_M)
        return _finite_result(self.nu @ self.flux_from_activities(c / self.standard_concentration_M), 'dc/dt')

    def chemical_potentials(self, activities):
        a = self._activities(activities, interior=True)
        return _finite_result(self.species_mu0 + self.rt * np.log(a), 'chemical potentials')

    def affinities(self, activities):
        """A=-nu.T mu in kcal/mol; restricted to positive activities."""
        return _finite_result(-(self.nu.T @ self.chemical_potentials(activities)), 'affinities')

    def ideal_free_energy_density(self, concentrations_M):
        """Sum c_i [mu0_i + RT(ln(c_i/c0)-1)], kcal/L; interior only."""
        c = self._concentrations(concentrations_M, interior=True)
        return float(_finite_result(np.dot(c, self.chemical_potentials(c / self.standard_concentration_M) - self.rt), 'free energy density'))

    def ideal_dissipation(self, concentrations_M):
        """Closed ideal isothermal interior diagnostic, not nonideal/chemostat proof."""
        c = self._concentrations(concentrations_M, interior=True)
        activities = c / self.standard_concentration_M
        flux = self.flux_from_activities(activities)
        affinity = self.affinities(activities)
        contributions = _finite_result(flux * affinity, 'dissipation')
        dGdt = float(np.dot(self.chemical_potentials(activities), self.rhs(c)))
        return dict(flux_M_s=flux, affinity_kcal_mol=affinity,
                    reaction_dissipation_kcal_L_s=contributions,
                    dG_dt_kcal_L_s=dGdt,
                    minus_flux_affinity_kcal_L_s=-float(contributions.sum()),
                    entropy_production_kcal_L_K_s=float(contributions.sum()) / self.temperature_K,
                    physical_rates_validated=False,
                    domain='ideal closed isothermal system; strictly positive concentrations')

    def shifted_energy_reference(self, element_shifts_kcal_mol: Mapping[str, float],
                                 charge_shift_kcal_mol=0.0):
        """Change conserved elemental/charge reference in species AND transition state."""
        if not set(element_shifts_kcal_mol) <= set(self._elements):
            raise ValueError('Unknown elemental reference shift')
        shifts = np.array([_scalar(element_shifts_kcal_mol.get(e, 0), 'element reference') for e in self._elements])
        delta = self._inventory.T @ shifts + self.charges * _scalar(charge_shift_kcal_mol, 'charge reference')
        return replace(self, species_mu0=self.species_mu0 + delta,
                       transition_mu0=self.transition_mu0 + self.alpha.T @ delta)


def observational_identifiability(observation_jacobian, parameter_scales,
                                  observation_noise_scales, relative_tolerance=1e-8):
    """Local scaled observation-Jacobian SVD, not global structural identifiability.

    Jacobian entries are d(observable)/d(parameter) at a supplied parameter point.
    Noise scales have observable units; parameter scales have parameter units.
    This function requires actual supplied derivatives; it does not infer them
    from the reaction network and does not establish identifiability without data.
    """
    raw = np.asarray(observation_jacobian, dtype=float)
    if raw.ndim != 2 or min(raw.shape) < 1:
        raise ValueError('A nonempty two-dimensional observation Jacobian is required')
    matrix = _array(raw, raw.shape, 'observation_jacobian')
    ps = _array(parameter_scales, (matrix.shape[1],), 'parameter_scales')
    noise = _array(observation_noise_scales, (matrix.shape[0],), 'observation_noise_scales')
    tolerance = _scalar(relative_tolerance, 'relative_tolerance', positive=True)
    if tolerance >= 1 or np.any(ps <= 0) or np.any(noise <= 0):
        raise ValueError('Positive scales and 0 < relative_tolerance < 1 required')
    scaled = _finite_result(matrix * ps[None, :] / noise[:, None], 'scaled Jacobian')
    _, singular, vt = np.linalg.svd(scaled, full_matrices=True)
    threshold = tolerance * (singular[0] if singular.size else 0)
    rank = int(np.sum(singular > threshold))
    condition = float(singular[0] / singular[-1]) if rank == matrix.shape[1] else None
    return dict(numerical_rank=rank, parameter_count=matrix.shape[1], observation_count=matrix.shape[0],
                singular_values=singular, threshold=threshold,
                locally_full_column_rank=rank == matrix.shape[1], condition_number=condition,
                scaled_null_directions=vt[rank:, :], global_structural_identifiability_established=False,
                physical_parameter_estimation_validated=False)
