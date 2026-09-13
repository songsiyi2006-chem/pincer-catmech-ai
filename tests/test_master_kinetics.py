"""Numerical and material-ledger verification; fixtures are not chemistry."""
import numpy as np
import pytest

from pincer_catmech.kinetics.master_kinetics import (
    MasterKinetics, MasterRates, MasterThermochemistry, S, METAL,
    BASE_EQUIVALENTS, SPECIES, STEP_NAMES, software_fixture_rates,
    software_fixture_initial,
)
from pincer_catmech.kinetics.microkinetics import EvidenceError


def model(temperature=383.15):
    return MasterKinetics(software_fixture_rates(temperature), tbuoh_activity=.1,
                          allow_software_fixture=True)


@pytest.mark.parametrize("seed", range(4))
def test_hand_jacobian_matches_central_difference(seed):
    m = model()
    c = np.random.default_rng(seed).uniform(.001, .8, 18)
    if seed == 0:
        c[:] = 0  # The polynomial derivative must work on the boundary too.
    h = 1e-6
    fd = np.column_stack([(m.rhs(0, c+h*np.eye(18)[j])-m.rhs(0, c-h*np.eye(18)[j]))/(2*h) for j in range(18)])
    assert np.max(np.abs(fd-m.jacobian(0, c))) < 1e-5
    np.testing.assert_allclose(m.rhs(0, c), S @ m.rates(c), atol=1e-13)
    np.testing.assert_allclose(METAL @ m.jacobian(0, c), 0, atol=1e-12)


def test_exact_stoichiometric_invariants_and_dimer_weight():
    assert len(SPECIES) == len(STEP_NAMES) == 18
    assert S.shape == (18, 18)
    np.testing.assert_array_equal(METAL @ S, 0)
    np.testing.assert_array_equal(BASE_EQUIVALENTS @ S, 0)
    assert METAL[6] == 2
    assert S[3, 13] == -2 and S[6, 13] == 1


@pytest.mark.parametrize("temperature,base", [(360., .01), (383.15, .05), (420., .2)])
def test_radau_500_points_conserves_inventory(temperature, base):
    out = model(temperature).integrate(software_fixture_initial(free_base_equiv=base), end_time=500.)
    assert out.concentrations.shape == (500, 18)
    assert out.metal_error_M < 1e-10
    assert out.base_equivalent_error_M < 1e-10
    assert out.minimum_concentration_M > -1e-10
    assert out.evidence_kind == "software_fixture"
    assert out.cleavage_fragment_formed_M[-1] >= 0
    assert out.active_metal_fraction[-1] < 1


def test_fixture_requires_explicit_opt_in_and_computed_needs_provenance():
    with pytest.raises(EvidenceError):
        MasterKinetics(software_fixture_rates(), tbuoh_activity=.1)
    with pytest.raises(EvidenceError):
        MasterRates(np.ones(18), np.zeros(18), 383.15, "computed")
    with pytest.raises(EvidenceError):
        MasterThermochemistry(383.15, {}, {}, "").validate()
    with pytest.raises(ValueError):
        MasterKinetics(software_fixture_rates(), tbuoh_activity=np.nan, allow_software_fixture=True)


def test_barrier_perturbation_preserves_forward_reverse_ratio():
    rates = software_fixture_rates()
    perturbed = rates.perturbed(4, .2)
    assert perturbed.forward[4]/perturbed.reverse[4] == pytest.approx(rates.forward[4]/rates.reverse[4])
    assert not rates.forward.flags.writeable
    with pytest.raises(ValueError):
        rates.perturbed(18, .2)


def test_tangent_control_matches_independent_barrier_perturbations():
    m = model()
    c0 = software_fixture_initial()
    trajectory = m.integrate(c0, end_time=100., points=500)
    control = m.control_coefficients(trajectory)
    assert control["rate_control"].shape == (500, 18)
    assert np.isnan(control["rate_control"][0]).all()
    assert control["steady_state_campbell_drc_established"] is False
    assert control["maximum_metal_sensitivity_residual_M"] < 1e-10
    # Distinct catalytic, shuttle and irreversible loss parameters.
    for step in (1, 10, 16):
        h = 1e-3
        results = []
        for theta in (h, -h):
            perturbed = MasterKinetics(m.parameters.perturbed(step, theta),
                                      tbuoh_activity=.1, allow_software_fixture=True)
            results.append(perturbed.integrate(c0, end_time=100., points=500))
        plus, minus = results
        fd = (np.log(plus.tof_per_second[-1])-np.log(minus.tof_per_second[-1]))/(2*h)
        assert control["rate_control"][-1, step] == pytest.approx(fd, abs=2e-5)
        sc = []
        for out in results:
            r = out.fluxes[-1]
            sc.append(np.log(r[5]/(r[5]+r[14]+r[15])))
        assert control["instantaneous_selectivity_control"][-1, step] == pytest.approx((sc[0]-sc[1])/(2*h), abs=2e-5)


def test_dimer_only_and_base_depletion_edge_cases():
    f, b = np.zeros(18), np.zeros(18)
    f[13], b[13] = 5000, .001
    m = MasterKinetics(MasterRates(f, b, 383.15, "software_fixture"), tbuoh_activity=0, allow_software_fixture=True)
    c0 = np.zeros(18)
    c0[3] = .01
    out = m.integrate(c0, end_time=10.)
    np.testing.assert_allclose(out.concentrations[:, 3]+2*out.concentrations[:, 6], .01, atol=1e-11)
    f[:] = 0
    f[16] = 100
    m = MasterKinetics(MasterRates(f, np.zeros(18), 383.15, "software_fixture"), tbuoh_activity=0, allow_software_fixture=True)
    c0 = np.zeros(18)
    c0[0], c0[17] = .01, .002
    out = m.integrate(c0, end_time=100.)
    assert out.concentrations[-1, 8] == pytest.approx(.002, abs=1e-10)
    assert out.concentrations[-1, 0] == pytest.approx(.008, abs=1e-10)
