"""Synthetic mathematical fixtures only: no physical catalyst rates or data."""
from dataclasses import replace
import numpy as np
import pytest

from pincer_catmech.kinetics.thermodynamic_network import (
    ThermodynamicNetwork, observational_identifiability)


def synthetic_cycle():
    return ThermodynamicNetwork(species=('A_fixture', 'B_fixture', 'C_fixture'),
        reactions=('AB_fixture', 'BC_fixture', 'CA_fixture'), species_mu0=[0., .3, -.2],
        transition_mu0=[21., 23., 24.], alpha=[[1,0,0],[0,1,0],[0,0,1]],
        beta=[[0,0,1],[1,0,0],[0,1,0]], compositions=({'H':1},)*3,
        charges=[0,0,0], temperature_K=383.15)


def synthetic_association(c0=1.):
    return ThermodynamicNetwork(species=('A_fixture', 'B_fixture', 'AB_fixture'),
        reactions=('association_fixture',), species_mu0=[-.1, .2, -.5],
        transition_mu0=[22.], alpha=[[1],[1],[0]], beta=[[0],[0],[1]],
        compositions=({'H':1},{'Cl':1},{'H':1,'Cl':1}),
        charges=[1,-1,0], temperature_K=383.15, standard_concentration_M=c0)


def test_cycle_detailed_balance_and_wegscheider():
    model=synthetic_cycle()
    lf,lr=model.log_frequency_constants()
    np.testing.assert_allclose(lf-lr,model.log_equilibrium_constants(),rtol=0,atol=5e-15)
    assert np.sum(model.log_equilibrium_constants()) == pytest.approx(0,abs=1e-15)
    assert np.prod(np.exp(model.log_equilibrium_constants())) == pytest.approx(1)
    activities=np.exp(-model.species_mu0/model.rt)
    np.testing.assert_allclose(model.flux_from_activities(activities),0,atol=1e-14)
    np.testing.assert_allclose(model.affinities(activities),0,atol=1e-15)
    assert model.physical_rates_validated is False


@pytest.mark.parametrize('seed',range(8))
def test_nonnegative_entropy_and_free_energy_derivative(seed):
    model=synthetic_association()
    c=np.random.default_rng(seed).uniform(.01,1.5,3)
    result=model.ideal_dissipation(c)
    assert np.all(result['reaction_dissipation_kcal_L_s']>=-1e-14)
    assert result['entropy_production_kcal_L_K_s']>=-1e-14
    assert result['dG_dt_kcal_L_s']<=1e-14
    assert result['dG_dt_kcal_L_s']==pytest.approx(result['minus_flux_affinity_kcal_L_s'],rel=1e-12,abs=1e-14)
    direction=model.rhs(c)
    dt=1e-6/max(1,np.linalg.norm(direction))
    difference=(model.ideal_free_energy_density(c+dt*direction)-model.ideal_free_energy_density(c-dt*direction))/(2*dt)
    assert difference==pytest.approx(result['dG_dt_kcal_L_s'],rel=1e-6,abs=1e-8)


@pytest.mark.parametrize('concentrations',[[0.,0.,0.],[0.,1.,1.],[1.,0.,1.],[1.,1.,0.],[0.,0.,1.]])
def test_boundary_positive_invariance_and_zero_power(concentrations):
    model=synthetic_association()
    c=np.array(concentrations); rhs=model.rhs(c)
    assert np.all(np.isfinite(rhs))
    assert np.all(rhs[c==0]>=0)
    assert np.dot(model.charges,rhs)==pytest.approx(0)
    np.testing.assert_allclose(np.array([[1,0,1],[0,1,1]])@rhs,0,atol=1e-15)
    if not np.any(c): np.testing.assert_array_equal(rhs,0)


def test_interior_diagnostics_reject_zero_without_clipping():
    model=synthetic_association()
    for function in [model.ideal_dissipation,model.ideal_free_energy_density,model.chemical_potentials,model.affinities]:
        with pytest.raises(ValueError,match='strictly positive'):
            function([1.,1.,0.])
    with pytest.raises(ValueError,match='nonnegative'):
        model.rhs([1.,1.,-1e-15])


def test_molecularity_and_standard_concentration_units():
    model=synthetic_association(c0=2.5)
    c=np.array([.3,.4,.1]); constants=model.concentration_rate_constants()
    assert constants['forward_molarity_power'][0]==-1
    assert constants['reverse_molarity_power'][0]==0
    expected=constants['forward'][0]*c[0]*c[1]-constants['reverse'][0]*c[2]
    assert model.flux_from_activities(c/2.5)[0]==pytest.approx(expected,rel=1e-13)
    # Changing c0 while preserving mu0 changes the specified thermodynamic model;
    # this check verifies units, not an illegitimate standard-state invariance.
    other=synthetic_association(c0=1.)
    k1=other.concentration_rate_constants()
    assert constants['forward'][0]==pytest.approx(k1['forward'][0]/2.5)
    assert constants['reverse'][0]==pytest.approx(k1['reverse'][0])


def test_element_and_charge_reference_gauge_invariance():
    model=synthetic_association(); shifted=model.shifted_energy_reference({'H':4.2,'Cl':-2.7},charge_shift_kcal_mol=.6)
    activities=np.array([.1,.3,.9])
    np.testing.assert_allclose(shifted.flux_from_activities(activities),model.flux_from_activities(activities),rtol=1e-12)
    np.testing.assert_allclose(shifted.affinities(activities),model.affinities(activities),atol=3e-15)
    np.testing.assert_allclose(shifted.log_equilibrium_constants(),model.log_equilibrium_constants(),atol=1e-15)
    assert shifted.ideal_dissipation(activities)['dG_dt_kcal_L_s']==pytest.approx(model.ideal_dissipation(activities)['dG_dt_kcal_L_s'])


def test_element_conservation_rejects_missing_atoms():
    with pytest.raises(ValueError,match='elemental conservation'):
        replace(synthetic_association(),compositions=({'H':1},{'Cl':1},{'H':1}))


def test_charge_conservation_rejects_hidden_counterion():
    with pytest.raises(ValueError,match='charge conservation'):
        replace(synthetic_association(),charges=[1,0,0])


@pytest.mark.parametrize('value',[-1,1.5,True])
def test_invalid_stoichiometry_rejected(value):
    with pytest.raises(ValueError):
        replace(synthetic_association(),alpha=[[value],[1],[0]])


@pytest.mark.parametrize('field,value',[
    ('temperature_K',0),('temperature_K',float('nan')),('temperature_K',True),
    ('standard_concentration_M',float('inf')),('standard_concentration_M',0),
    ('species_mu0',[0,float('nan'),0]),('transition_mu0',[float('inf')]),
    ('charges',[0,.5,0]),('compositions',({'Fake':1},{'Cl':1},{'H':1,'Cl':1})),
    ('species',('A','A','C')),('beta',[[1],[1],[0]])])
def test_invalid_model_inputs_rejected(field,value):
    with pytest.raises(ValueError):
        replace(synthetic_association(),**{field:value})


@pytest.mark.parametrize('values',[[1,np.nan,1],[1,np.inf,1],[1,-.1,1]])
def test_nonfinite_or_negative_states_rejected(values):
    with pytest.raises(ValueError): synthetic_association().flux_from_activities(values)


def test_immutable_input_copy():
    model=synthetic_association()
    with pytest.raises(ValueError): model.alpha.setflags(write=True)
    with pytest.raises(TypeError): model.compositions[0]['H']=2


def test_local_rank_only_detects_collinear_observations():
    result=observational_identifiability([[1,2],[2,4],[3,6]],[1,1],[1,1,1])
    assert result['numerical_rank']==1
    assert result['scaled_null_directions'].shape==(1,2)
    assert not result['locally_full_column_rank']
    assert result['condition_number'] is None
    assert not result['global_structural_identifiability_established']
    assert not result['physical_parameter_estimation_validated']


def test_rank_is_invariant_to_consistent_observation_and_parameter_units():
    raw=np.array([[1.,0.],[0.,2.],[1.,1.]])
    first=observational_identifiability(raw,[2,3],[.1,.2,.3])
    converted=observational_identifiability(raw*1000/60,[2*60,3*60],np.array([.1,.2,.3])*1000)
    np.testing.assert_allclose(first['singular_values'],converted['singular_values'])
    assert first['numerical_rank']==converted['numerical_rank']==2


def test_underdetermined_nullspace_includes_unobserved_dimensions():
    result=observational_identifiability([[1.,2.,3.]],[1,1,1],[1])
    assert result['numerical_rank']==1
    assert result['scaled_null_directions'].shape==(2,3)


@pytest.mark.parametrize('jac,ps,noise',[
    ([[np.nan]],[1],[1]),([[1]],[0],[1]),([[1]],[1],[0]),([],[],[]),([[1]],[np.inf],[1])])
def test_invalid_identifiability_inputs(jac,ps,noise):
    with pytest.raises(ValueError): observational_identifiability(jac,ps,noise)
