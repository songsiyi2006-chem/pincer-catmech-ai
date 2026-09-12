"""Analytical/software cases only: no fake computational campaign observations."""
from dataclasses import replace
import hashlib
import math

import numpy as np
import pytest

from pincer_catmech.kinetics.microkinetics import (
    BALANCES, SPECIES, ComputedFreeEnergy, EvidenceError, MicrokineticModel,
    ReactionChannel, ThermodynamicSnapshot, STEPS, build_model, campbell_drc, campaign_readiness,
    eyring_rate_constant, initial_concentrations, run_condition_grid, solve_batch,
)
from pincer_catmech.models.surrogate_regressor import (
    BarrierObservation, CandidatePerformance, pareto_front, save_surrogate,
    train_surrogate, training_readiness,
)


def analytical_model(k=0.2, reverse=0.0):
    return MicrokineticModel((ReactionChannel("analytic_release", ((9, 1),), ((10, 1), (1, 1)), k, reverse),),
                            400.0, "analytical_unit_test")


def analytical_initial():
    state = dict.fromkeys(SPECIES, 0.0)
    state["complex_2"] = 0.01
    return state


@pytest.mark.parametrize("method", ["Radau", "BDF"])
def test_eight_ode_solution_against_exact_first_order_release(method):
    result = solve_batch(analytical_model(), analytical_initial(), end_time_seconds=10, method=method)
    expected = 0.01*(1-np.exp(-0.2*result.time_seconds))
    np.testing.assert_allclose(result.concentrations_mol_l[:, 10], expected, atol=2e-9, rtol=1e-6)
    assert len(result.integrated_species) == 8
    assert result.concentrations_mol_l.shape == (101, 11)
    assert result.conservation_max_abs_mol_l < 1e-14
    np.testing.assert_allclose(result.concentrations_mol_l @ BALANCES.T, np.full((101, 3), 0.01), atol=1e-14)
    assert result.average_tof_per_second == pytest.approx((1-math.exp(-2))/10, rel=1e-6)
    assert result.instantaneous_tof_per_second == pytest.approx(0.2*math.exp(-2), rel=1e-6)
    assert result.evidence_level == "analytical_unit_test"


def test_eyring_molecularity_and_reversible_detailed_balance():
    temperature, c0, forward, reaction = 400.0, 0.5, 16.0, -4.0
    k_forward = eyring_rate_constant(forward, temperature, 2, c0)
    k_reverse = eyring_rate_constant(forward-reaction, temperature, 1, c0)
    gas_constant = 8.31446261815324/4184
    assert k_forward/k_reverse == pytest.approx(math.exp(-reaction/(gas_constant*temperature))/c0)
    assert eyring_rate_constant(16, temperature, 2, c0) == pytest.approx(2*eyring_rate_constant(16, temperature, 1, c0))


@pytest.mark.parametrize("barrier,temperature,molecularity", [(-1, 400, 1), (12, 0, 1), (12, 400, 0), (float("nan"), 400, 1)])
def test_eyring_rejects_invalid_physical_input(barrier, temperature, molecularity):
    with pytest.raises(ValueError):
        eyring_rate_constant(barrier, temperature, molecularity)


def test_finite_time_drc_matches_analytical_sensitivity():
    result = campbell_drc(analytical_model(), analytical_initial(), end_time_seconds=2,
                          perturbation_kcal_mol=0.001)
    expected = 0.4/math.expm1(0.4)
    assert result[0].value == pytest.approx(expected, abs=2e-6)
    assert result[0].definition == "finite_time_apparent_DRC_of_average_TOF"


def test_zero_flux_drc_is_undefined():
    result = campbell_drc(analytical_model(k=0), analytical_initial(), end_time_seconds=2)
    assert result[0].value is None
    assert result[0].status == "undefined_nonpositive_net_TOF"


def test_invalid_network_cannot_destroy_catalyst_or_fragments():
    with pytest.raises(ValueError, match="conservation"):
        MicrokineticModel((ReactionChannel("bad", ((9, 1),), ((10, 1),), 0.1, 0),), 400, "analytical_unit_test")


def test_base_spectator_scales_both_forward_and_backward_fluxes():
    channel = ReactionChannel("base_path", ((6, 1),), ((7, 1), (8, 1)), 3.0, 2.0, 0.5)
    concentrations = np.ones(11)
    assert channel.flux(concentrations) == pytest.approx(0.5)
    assert replace(channel, base_activity=2).flux(concentrations) == pytest.approx(2)
    assert channel.k_forward/channel.k_reverse == replace(channel, base_activity=2).k_forward/replace(channel, base_activity=2).k_reverse


def test_empty_computed_campaign_is_gated_not_filled():
    readiness = campaign_readiness({})
    assert readiness["status"] == "insufficient_computed_data"
    assert len(readiness["reasons"]) == 11
    with pytest.raises(EvidenceError, match="Missing"):
        run_condition_grid({}, end_time_seconds=60)
    snapshot = ThermodynamicSnapshot(400, {}, {})
    with pytest.raises(EvidenceError, match="eleven"):
        snapshot.validate()


def test_computed_record_rejects_test_data_even_with_matching_hash(tmp_path):
    source = tmp_path / "analytical-fixture.txt"
    source.write_text("Explicit analytical unit-test fixture, not quantum output", encoding="utf-8")
    record = ComputedFreeEnergy(1, 400, str(source), hashlib.sha256(source.read_bytes()).hexdigest(),
                                "analytical", True, True, evidence_kind="analytical_unit_test")
    with pytest.raises(EvidenceError, match="computed"):
        record.validate(400)


def test_missing_hash_and_unaccepted_ts_are_rejected(tmp_path):
    record = ComputedFreeEnergy(12, 400, str(tmp_path/"absent.out"), "0"*64,
                                "GFN2-xTB", True, True, (-500,), True)
    with pytest.raises(EvidenceError, match="artifact"):
        record.validate(400, transition_state=True)
    with pytest.raises(EvidenceError, match="accepted"):
        replace(record, accepted=False).validate(400, transition_state=True)


def test_empty_surrogate_data_returns_insufficient_data_without_sklearn(tmp_path):
    assert training_readiness([])["status"] == "insufficient_computed_data"
    outcome = train_surrogate([])
    assert outcome.model is None
    assert outcome.status == "insufficient_computed_data"
    with pytest.raises(EvidenceError, match="No trained"):
        save_surrogate(outcome, tmp_path)
    assert not list(tmp_path.iterdir())


def test_pareto_refuses_analytical_tof_and_inconsistent_conditions():
    kinetics = solve_batch(analytical_model(), analytical_initial(), end_time_seconds=1)
    row = CandidatePerformance("unit", "Mn", kinetics, 400, 1.0, 1.0)
    with pytest.raises(EvidenceError, match="accepted"):
        pareto_front([row])
    with pytest.raises(ValueError, match="share"):
        pareto_front([row, replace(row, catalyst_id="unit2", temperature=410)])
    assert pareto_front([]) == []


def test_grouped_cv_algorithm_on_explicit_analytical_unit_cases(monkeypatch):
    pytest.importorskip("sklearn")
    # Bypass scientific record validation only inside this numerical algorithm
    # test. No artifact is saved and these labels never enter campaign datasets.
    monkeypatch.setattr(BarrierObservation, "validate", lambda self, names: self.descriptors["x"]*2+1)
    rows = []
    for group in range(5):
        for metal in range(2):
            record = ComputedFreeEnergy(0, 400, "analytical_only", f"{group}-{metal}",
                                        "analytical_test", False, False, evidence_kind="analytical_unit_test")
            rows.append(BarrierObservation(f"test-{group}-{metal}", f"ligand{group}", "R",
                                           {"x": group+metal/10}, record, ((record, 1),)))
    outcome = train_surrogate(rows, feature_names=("x",), ensemble_members=3)
    assert outcome.model is not None
    assert len(outcome.fold_records) == 5
    for fold in outcome.fold_records:
        assert not set(fold["train_ligand_groups"]) & set(fold["test_ligand_groups"])
        assert not set(fold["train_catalysts"]) & set(fold["test_catalysts"])
    prediction = outcome.model.predict([{"x": 2}, {"x": 100}])
    assert prediction["outside_training_box"] == [False, True]
    assert all(value >= 0 for value in prediction["ensemble_std_kcal_mol"])
    assert "not_calibrated" in prediction["uncertainty_kind"]


def test_loadings_are_mol_percent_relative_to_alcohol():
    initial = initial_concentrations(alcohol_mol_l=0.5, catalyst_loading_mol_percent=2)
    assert initial["cat"] == pytest.approx(0.01)


def _algorithm_snapshot(monkeypatch, *, assisted=True, co_shuttle=True):
    """Explicit analytical algebra fixture; never an accepted campaign dataset.

    Only the file/quantum-review validator is bypassed here. Composition,
    protocol, missing-step gates, Eyring construction and integration execute.
    """
    from collections import Counter
    monkeypatch.setattr(ComputedFreeEnergy, "validate", lambda *args, **kwargs: None)

    def combine(*parts):
        total = Counter()
        for part in parts:
            total.update(part)
        return dict(total)

    catalyst = {"Ru": 1, "C": 10, "H": 10, "N": 1, "P": 2, "O": 1}
    alcohol, amine = {"C": 7, "H": 8, "O": 1}, {"C": 6, "H": 7, "N": 1}
    aldehyde, water = {"C": 7, "H": 6, "O": 1}, {"H": 2, "O": 1}
    imine, product = {"C": 13, "H": 11, "N": 1}, {"C": 13, "H": 13, "N": 1}
    hydrogenated = combine(catalyst, {"H": 2})
    compositions = {"cat": catalyst, "alcohol": alcohol, "amine": amine,
        "complex_1": combine(catalyst, alcohol), "aldehyde": aldehyde,
        "cat_h2": hydrogenated, "hemiaminal": combine(aldehyde, amine),
        "imine": imine, "water": water, "complex_2": combine(imine, hydrogenated), "product": product}

    def record(name, composition, energy=0, charge=0, ts=False):
        return ComputedFreeEnergy(energy, 400, "analytical_unit_test_only", name,
            "analytical_algorithm_test", False, False, (-500,) if ts else (), False,
            evidence_kind="analytical_unit_test", solvent="test_medium", solvation_state="test_reference",
            composition=tuple(composition.items()), charge=charge)

    species = {name: record(name, composition) for name, composition in compositions.items()}
    states = {step: record(step, combine(*(compositions[name] for name in reactants)), 20, ts=True)
              for step, (reactants, _) in STEPS.items()}
    base = record("base", {"C": 4, "H": 9, "O": 1}, -10, -1) if assisted else None
    conjugate = record("conjugate", {"C": 4, "H": 10, "O": 1}, -5) if assisted and co_shuttle else None
    shuttle_parts = [dict(base.composition)] if assisted else []
    if conjugate:
        shuttle_parts.append(dict(conjugate.composition))
    paths = {"dehydration": record("dehydration_base", combine(compositions["hemiaminal"], *shuttle_parts),
              20 + base.gibbs_kcal_mol + (conjugate.gibbs_kcal_mol if conjugate else 0), -1, True)} if assisted else None
    return ThermodynamicSnapshot(400, species, states, paths, base, conjugate)


def test_all_six_steps_and_two_species_shuttle_keep_composition_and_detailed_balance(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch)
    model = build_model(snapshot, base_concentration_mol_l=0.05, co_shuttle_concentration_mol_l=0.1)
    assert len(model.channels) == 7
    assert np.array_equal(BALANCES @ model.stoichiometry, np.zeros((3, 7)))
    base = next(channel for channel in model.channels if channel.name.endswith(":base"))
    normal = next(channel for channel in model.channels if channel.name == "dehydration")
    assert base.base_activity == pytest.approx(0.005)
    assert base.k_forward == pytest.approx(normal.k_forward)
    assert base.k_reverse == pytest.approx(normal.k_reverse)
    concentration = np.ones(11)
    concentration[6] = 2
    assert base.flux(concentration) == pytest.approx(0.005*normal.flux(concentration))
    for channel in model.channels:
        assert channel.k_forward/channel.k_reverse == pytest.approx(1)


def test_model_rates_have_zero_net_flux_at_thermodynamic_equilibrium(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch, assisted=False)
    energies = np.linspace(-0.4, 0.6, 11)
    snapshot = replace(snapshot, species={name: replace(snapshot.species[name], gibbs_kcal_mol=float(energy))
                                        for name, energy in zip(SPECIES, energies)})
    model = build_model(snapshot)
    rt = (8.31446261815324/4184)*snapshot.temperature
    equilibrium = np.exp(-energies/rt)  # c0=1 M and zero conserved-reference potentials.
    for channel in model.channels:
        reaction_g = sum(n*energies[i] for i, n in channel.products)-sum(n*energies[i] for i, n in channel.reactants)
        assert channel.k_forward/channel.k_reverse == pytest.approx(math.exp(-reaction_g/rt))
        assert channel.flux(equilibrium) == pytest.approx(0, abs=1e-12)


def test_computed_channel_cannot_drop_water_while_preserving_fragment_balances(monkeypatch):
    model = build_model(_algorithm_snapshot(monkeypatch, assisted=False))
    channels = tuple(replace(channel, products=((7, 1),)) if channel.name == "dehydration" else channel
                     for channel in model.channels)
    with pytest.raises(EvidenceError, match="stoichiometry"):
        replace(model, channels=channels)


@pytest.mark.parametrize("failure", ["missing_step", "mixed_solvent", "wrong_charge", "missing_shuttle_atom"])
def test_snapshot_rejects_missing_or_inconsistent_reaction_evidence(monkeypatch, failure):
    snapshot = _algorithm_snapshot(monkeypatch)
    if failure == "missing_step":
        snapshot = replace(snapshot, transition_states={name: value for name, value in snapshot.transition_states.items()
                                                       if name != "hydrogenation"})
    elif failure == "mixed_solvent":
        snapshot = replace(snapshot, species={**snapshot.species, "water": replace(snapshot.species["water"], solvent="gas")})
    else:
        state = snapshot.base_assisted_transition_states["dehydration"]
        state = replace(state, charge=0) if failure == "wrong_charge" else replace(state, composition=(("C", 1),))
        snapshot = replace(snapshot, base_assisted_transition_states={"dehydration": state})
    with pytest.raises(EvidenceError):
        build_model(snapshot, base_concentration_mol_l=0.05, co_shuttle_concentration_mol_l=0.1)


def test_base_inputs_require_explicit_path_and_shuttle_reference(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch, assisted=False)
    with pytest.raises(EvidenceError, match="base-assisted"):
        build_model(snapshot, base_concentration_mol_l=0.05)
    with pytest.raises(EvidenceError, match="co-shuttle"):
        build_model(snapshot, co_shuttle_concentration_mol_l=0.1)


def test_condition_grid_records_total_and_free_base_separately(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch)
    arguments = dict(end_time_seconds=0.1, temperatures=(400,), loadings_mol_percent=(1,),
        base_equivalents=(0.05,), alcohol_mol_l=0.5, calculate_drc=False, co_shuttle_concentration_mol_l=0.1)
    with pytest.raises(EvidenceError, match="free-base fraction"):
        run_condition_grid({400: snapshot}, **arguments)
    row, = run_condition_grid({400: snapshot}, free_base_fraction=0.2, **arguments)
    assert row["total_base_concentration_mol_l"] == pytest.approx(0.025)
    assert row["free_base_concentration_mol_l"] == pytest.approx(0.005)
    assert row["result"].free_base_concentration_mol_l == pytest.approx(0.005)
    assert row["result"].co_shuttle_concentration_mol_l == 0.1
    assert row["result"].initial_concentrations_mol_l[0] == 0.5
    assert row["result"].temperature == 400


def test_declaring_computed_model_without_snapshot_cannot_bypass_gate():
    with pytest.raises(EvidenceError, match="source snapshot"):
        replace(analytical_model(), evidence_level="computed_accepted_thermochemistry")


def test_pareto_rejects_different_free_base_or_incorrect_declared_temperature(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch)
    model = build_model(snapshot, base_concentration_mol_l=0.01, co_shuttle_concentration_mol_l=0.1)
    kinetics = solve_batch(model, initial_concentrations(), end_time_seconds=0.1)
    row = CandidatePerformance("algorithm-only-1", "Mn", kinetics, 400, 1, 0.05)
    with pytest.raises(ValueError, match="share"):
        pareto_front([row, replace(row, catalyst_id="algorithm-only-2",
                                 kinetics=replace(kinetics, free_base_concentration_mol_l=0.02))])
    with pytest.raises(EvidenceError, match="temperature"):
        pareto_front([replace(row, temperature=410)])
    with pytest.raises(EvidenceError, match="loading"):
        pareto_front([replace(row, loading_mol_percent=2)])


def test_matching_artifact_hash_does_not_replace_identity_or_solvent_review(tmp_path):
    source = tmp_path/"deliberately_incomplete_evidence.txt"
    source.write_text("Unit rejection case; no quantum chemistry was performed.")
    # These claimed flags are deliberately insufficient and must be rejected.
    record = ComputedFreeEnergy(1, 400, str(source), hashlib.sha256(source.read_bytes()).hexdigest(),
                                "claimed_method", True, True)
    with pytest.raises(EvidenceError, match="solvent"):
        record.validate(400)
    with pytest.raises(EvidenceError, match="identity"):
        replace(record, solvent="toluene", solvation_state="gsolv").validate(400)
    with pytest.raises(EvidenceError, match="reaction motion"):
        replace(record, solvent="toluene", solvation_state="gsolv", chemical_identity_verified=True,
                composition=(("H", 2),), charge=0, imaginary_frequencies_cm1=(-500,),
                connectivity_verified=True).validate(400, transition_state=True)


def test_surrogate_rejects_an_incomplete_stoichiometric_reactant_reference(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch, assisted=False)
    row = BarrierObservation("analytical-case", "test-ligand", "R", {"x": 1},
                             snapshot.transition_states["dehydrogenation"], ((snapshot.species["water"], 1),))
    with pytest.raises(EvidenceError, match="elemental composition"):
        row.validate(("x",))


def test_surrogate_does_not_mix_bound_complex_and_separated_reference_labels(monkeypatch):
    snapshot = _algorithm_snapshot(monkeypatch, assisted=False)
    bound = BarrierObservation("bound-reference-test", "ligand1", "R", {"x": 1},
        snapshot.transition_states["dehydrogenation"], ((snapshot.species["complex_1"], 1),))
    separated = replace(bound, catalyst_id="separated-reference-test", backbone="ligand2",
        transition_state=replace(bound.transition_state, source_sha256="distinct-test-label"),
        reactants=((snapshot.species["cat"], 1), (snapshot.species["alcohol"], 1)),
        reference_convention="separated_reactants")
    assert bound.validate(("x",)) == 20
    assert separated.validate(("x",)) == 20
    readiness = training_readiness([bound, separated], feature_names=("x",))
    assert any("activation-reference" in reason for reason in readiness["reasons"])
    with pytest.raises(EvidenceError, match="molecularity"):
        replace(separated, reference_convention="preassociated_complex").validate(("x",))


def test_drc_perturbs_forward_and_reverse_rates_by_the_same_factor(monkeypatch):
    from pincer_catmech.kinetics import microkinetics as module
    original_solver = module.solve_batch
    observed = []

    def observed_solver(model, *args, **kwargs):
        observed.append((model.channels[0].k_forward, model.channels[0].k_reverse))
        return original_solver(model, *args, **kwargs)

    monkeypatch.setattr(module, "solve_batch", observed_solver)
    campbell_drc(analytical_model(k=0.2, reverse=0.01), analytical_initial(), end_time_seconds=2)
    assert len(observed) == 2
    assert all(forward/reverse == pytest.approx(20) for forward, reverse in observed)
