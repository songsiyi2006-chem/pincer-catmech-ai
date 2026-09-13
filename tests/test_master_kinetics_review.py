"""Independent adversarial tests. Algebra-only records never leave this file.

The test-only monkeypatch explicitly disables quantum-record attestation when
exercising mathematical construction. It does not produce computed chemistry.
"""
from collections import Counter
from dataclasses import replace
import hashlib

import numpy as np
import pytest

from pincer_catmech.kinetics.master_kinetics import (
    BASE_EQUIVALENTS, INDEX, METAL, REACTIONS, S, SPECIES, STEP_NAMES,
    MasterKinetics, MasterRates, MasterThermochemistry,
    software_fixture_initial, software_fixture_rates,
)
from pincer_catmech.kinetics.microkinetics import ComputedFreeEnergy, EvidenceError, R_KCAL


def _algebra_snapshot(monkeypatch):
    """Known integer ledger, with no pretence that these are QM calculations."""
    monkeypatch.setattr(ComputedFreeEnergy, "validate", lambda *args, **kwargs: None)

    def add(*parts):
        result = Counter()
        for part in parts:
            result.update(part)
        return dict(result)

    cat = dict(Ru=1, C=10, H=20, N=1, P=2, O=1)
    alcohol, aldehyde = dict(C=7, H=8, O=1), dict(C=7, H=6, O=1)
    amine, imine, product = dict(C=6, H=7, N=1), dict(C=13, H=11, N=1), dict(C=13, H=13, N=1)
    hydrogenated = add(cat, dict(H=2))
    formulas = {
        "cat": cat, "cat_alcohol": add(cat, alcohol),
        "cat_h2_aldehyde": add(hydrogenated, aldehyde), "cat_h2": hydrogenated,
        "cat_h2_imine": add(hydrogenated, imine), "cat_product": add(cat, product),
        "dimer": add(hydrogenated, hydrogenated), "cat_co": add(cat, dict(C=1, O=1)),
        "cat_degraded": dict(Ru=1, C=9, H=17, N=1, P=2, O=1),
        "alcohol": alcohol, "aldehyde": aldehyde, "amine": amine,
        "hemiaminal": add(aldehyde, amine), "imine": imine, "product": product,
        "water": dict(H=2, O=1), "benzene": dict(C=6, H=6),
        "free_base": dict(C=4, H=9, O=1), "tbuoh": dict(C=4, H=10, O=1),
        "cleavage_fragment": dict(C=5, H=12, O=1),
    }

    def record(name, formula, *, charge=0, energy=0., ts=False):
        return ComputedFreeEnergy(
            energy, 400., "analytical_algebra_ONLY", name, "analytical_algebra_ONLY",
            False, False, (-500.,) if ts else (), False,
            evidence_kind="analytical_unit_test", solvent="test", solvation_state="test",
            composition=tuple(formula.items()), charge=charge,
        )

    states = {name: record(name, formula, charge=-1 if name in {"free_base", "cat_degraded"} else 0)
              for name, formula in formulas.items()}
    transitions = {}
    for j, (name, (left, _)) in enumerate(zip(STEP_NAMES, REACTIONS)):
        total = Counter()
        charge = 0
        for species, n in left.items():
            total.update({symbol: count*n for symbol, count in formulas[species].items()})
            charge += n*states[species].charge
        transitions[name] = record(name, dict(total), charge=charge, energy=24.+j*.01, ts=True)
    return MasterThermochemistry(400., states, transitions, "Algebra-only irreversible-sink unit test")


def test_computed_rate_constructor_rejects_arbitrary_arrays_and_wrong_temperature(monkeypatch):
    snapshot = _algebra_snapshot(monkeypatch)
    rates = snapshot.rates()
    with pytest.raises(EvidenceError, match="Eyring"):
        replace(rates, forward=rates.forward*2)
    with pytest.raises(EvidenceError, match="temperatures"):
        replace(rates, temperature=410.)
    with pytest.raises(EvidenceError, match="explicit evidence label"):
        replace(rates, ts_perturbations_theta=(.1,)+(0.,)*17)


def test_computed_ts_perturbation_has_explicit_label_and_source_binding(monkeypatch):
    rates = _algebra_snapshot(monkeypatch).rates()
    shifted = rates.perturbed(2, .3)
    assert shifted.evidence_kind == "computed_TS_perturbation"
    assert shifted.thermochemistry is rates.thermochemistry
    assert shifted.ts_perturbations_theta[2] == .3
    assert shifted.forward[2]/shifted.reverse[2] == pytest.approx(rates.forward[2]/rates.reverse[2])
    with pytest.raises(EvidenceError, match="Eyring"):
        replace(shifted, reverse=shifted.reverse*.8)
    restored = shifted.perturbed(2, -.3)
    assert restored.evidence_kind == "computed"
    np.testing.assert_allclose(restored.forward, rates.forward, rtol=1e-14)


def test_thermochemistry_mapping_cannot_change_after_validation(monkeypatch):
    snapshot = _algebra_snapshot(monkeypatch)
    with pytest.raises(TypeError):
        snapshot.states["cat"] = snapshot.states["water"]
    with pytest.raises(TypeError):
        snapshot.transition_states["alcohol_h_transfer"] = snapshot.states["cat"]


@pytest.mark.parametrize("step", STEP_NAMES)
def test_every_transition_state_is_required(monkeypatch, step):
    snapshot = _algebra_snapshot(monkeypatch)
    states = dict(snapshot.transition_states)
    del states[step]
    with pytest.raises(EvidenceError, match="all 18"):
        replace(snapshot, transition_states=states).rates()


def test_all_parallel_condensation_paths_share_thermodynamic_ratio(monkeypatch):
    snapshot = _algebra_snapshot(monkeypatch)
    states = dict(snapshot.states)
    states["water"] = replace(states["water"], gibbs_kcal_mol=1.2)
    rates = replace(snapshot, states=states).rates()
    expected = np.exp(-1.2/(R_KCAL*400.))
    for j in (7, 9, 10, 11, 12):
        assert rates.forward[j]/rates.reverse[j] == pytest.approx(expected)
    assert rates.forward[6]/rates.reverse[6] == pytest.approx(rates.forward[8]/rates.reverse[8])


def test_dimer_four_hydrogens_and_eliminated_fragment_close_atom_charge_ledger(monkeypatch):
    snapshot = _algebra_snapshot(monkeypatch)
    snapshot.validate()
    cat_h = dict(snapshot.states["cat"].composition)["H"]
    assert dict(snapshot.states["dimer"].composition)["H"] == 2*cat_h+4
    labels, ledger = snapshot.elemental_charge_ledger()
    np.testing.assert_array_equal(ledger @ S, 0)
    assert ledger[labels.index("charge"), INDEX["cat_degraded"]] == -1
    assert ledger[labels.index("charge"), INDEX["free_base"]] == -1
    # A partially degraded initial catalyst does not imply a newly formed fragment.
    initial = software_fixture_initial()
    initial[0] -= .001
    initial[8] += .001
    model = MasterKinetics(snapshot.rates(), tbuoh_activity=.1)
    trajectory = model.integrate(initial, end_time=5., points=11)
    assert trajectory.cleavage_fragment_formed_M[0] == 0
    assert trajectory.elemental_charge_error_M < 1e-10
    np.testing.assert_allclose(ledger @ (trajectory.concentrations-initial).T, 0, atol=1e-10)


@pytest.mark.parametrize("name", ["cat_h2", "dimer", "cat_degraded", "free_base", "cleavage_fragment", "tbuoh"])
def test_invalid_atom_or_charge_ledger_is_rejected(monkeypatch, name):
    snapshot = _algebra_snapshot(monkeypatch)
    states = dict(snapshot.states)
    record = states[name]
    states[name] = replace(record, charge=record.charge+1)
    with pytest.raises(EvidenceError):
        replace(snapshot, states=states).rates()


def test_single_channel_dimer_matches_closed_form_not_another_ode_solver():
    forward = np.zeros(18)
    forward[13] = 12.5
    model = MasterKinetics(MasterRates(forward, np.zeros(18), 400., "software_fixture"),
                          tbuoh_activity=0, allow_software_fixture=True)
    initial = np.zeros(18)
    initial[3] = .01
    trajectory = model.integrate(initial, end_time=30., points=41)
    expected = initial[3]/(1+2*forward[13]*initial[3]*trajectory.time)
    np.testing.assert_allclose(trajectory.concentrations[:, 3], expected, rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(trajectory.concentrations[:, 6], (initial[3]-expected)/2, rtol=1e-8, atol=1e-12)


def test_tangent_rejects_trajectory_from_different_rates_or_reservoir():
    rates = software_fixture_rates()
    model = MasterKinetics(rates, tbuoh_activity=.1, allow_software_fixture=True)
    trajectory = model.integrate(software_fixture_initial(), end_time=1., points=2)
    for other_rates, activity in ((rates.perturbed(0, .1), .1), (rates, .2)):
        other = MasterKinetics(other_rates, tbuoh_activity=activity, allow_software_fixture=True)
        with pytest.raises(ValueError, match="same rates"):
            other.control_coefficients(trajectory)


@pytest.mark.parametrize("temperature", [0., -1., np.nan, np.inf])
def test_fixture_rejects_nonphysical_temperature(temperature):
    with pytest.raises(ValueError, match="temperature"):
        software_fixture_rates(temperature)


@pytest.mark.parametrize("value", [0., -1., np.nan, np.inf])
def test_solver_rejects_invalid_tolerances_before_integration(value):
    model = MasterKinetics(software_fixture_rates(), tbuoh_activity=.1, allow_software_fixture=True)
    with pytest.raises(ValueError, match="tolerances"):
        model.integrate(software_fixture_initial(), rtol=value, end_time=1.)


def _attestation_shape(tmp_path):
    """Shape-valid attestation metadata, used only for rejection tests.

    A matching hash certifies artifact integrity, not that the content contains
    the asserted QM result. This fixture is never accepted for a calculation.
    """
    path = tmp_path / "NOT_A_CALCULATION.txt"
    path.write_text("Adversarial metadata test only; no quantum chemistry was run.", encoding="utf-8")
    return ComputedFreeEnergy(
        24., 400., str(path), hashlib.sha256(path.read_bytes()).hexdigest(),
        "test_attestation_shape", True, True, (-500.,), True,
        solvent="gas", solvation_state="none", chemical_identity_verified=True,
        reaction_mode_verified=True, composition=(("C", 1), ("H", 4)), charge=0,
    )


@pytest.mark.parametrize("changes", [
    {"evidence_kind": "software_fixture"}, {"accepted": False},
    {"stationarity_verified": False}, {"chemical_identity_verified": False},
    {"reaction_mode_verified": False}, {"connectivity_verified": False},
    {"source_sha256": "0"*64}, {"gibbs_kcal_mol": np.nan},
    {"temperature": 401.}, {"standard_state": "1atm"},
    {"imaginary_frequencies_cm1": ()}, {"imaginary_frequencies_cm1": (-500., -30.)},
    {"imaginary_frequencies_cm1": (500.,)}, {"imaginary_frequencies_cm1": (-50.,)},
    {"composition": (("C", True),)}, {"charge": True},
])
def test_computed_free_energy_adversarial_rejections(tmp_path, changes):
    record = replace(_attestation_shape(tmp_path), **changes)
    with pytest.raises(EvidenceError):
        record.validate(400., transition_state=True, dehydrogenation=True)


def test_source_modification_after_attestation_is_rejected(tmp_path):
    record = _attestation_shape(tmp_path)
    from pathlib import Path
    Path(record.source_path).write_text("Changed bytes", encoding="utf-8")
    with pytest.raises(EvidenceError, match="SHA256"):
        record.validate(400., transition_state=True)
