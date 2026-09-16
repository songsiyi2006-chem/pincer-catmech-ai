"""Synthetic mathematical audit; all species/energies are fixtures, not catalyst data."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'src'))
from pincer_catmech.kinetics.thermodynamic_network import ThermodynamicNetwork, observational_identifiability


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=REPO/'data/phase4/thermodynamic_fixture_audit_001')
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Fresh audit destination required')
    alpha=np.eye(3,dtype=int)
    beta=np.array([[0,0,1],[1,0,0],[0,1,0]])
    network=ThermodynamicNetwork(species=('A_fixture','B_fixture','C_fixture'),
        reactions=('A_to_B_fixture','B_to_C_fixture','C_to_A_fixture'),
        species_mu0=np.array([0.,1.,-1.]),transition_mu0=np.array([22.,23.,22.5]),
        alpha=alpha,beta=beta,compositions=({'C':1},{'C':1},{'C':1}),
        charges=np.array([0,0,0]),temperature_K=383.15)
    rng=np.random.default_rng(20260916)
    states=10**rng.uniform(-6,0,size=(200,3))
    defects=[]
    entropy=[]
    mass=[]
    gauge=[]
    shifted=network.shifted_energy_reference({'C':7.5})
    for c in states:
        d=network.ideal_dissipation(c)
        defects.append(abs(d['dG_dt_kcal_L_s']-d['minus_flux_affinity_kcal_L_s']))
        entropy.append(d['entropy_production_kcal_L_K_s'])
        mass.append(abs(network.rhs(c).sum()))
        gauge.append(float(np.max(np.abs(network.rhs(c)-shifted.rhs(c)))))
    weights=np.exp(-network.species_mu0/network.rt)
    eq=0.17*weights/weights.sum()
    equilibrium_flux=float(np.max(np.abs(network.flux_from_activities(eq))))
    cycle_logK_sum=float(np.sum(network.log_equilibrium_constants()))
    inward=[]
    for i in range(3):
        c=np.array([.1,.05,.02]);c[i]=0
        inward.append(float(network.rhs(c)[i]))
    rank1=observational_identifiability(np.array([[1.,1.],[2.,2.],[3.,3.]]),np.ones(2),np.ones(3))
    rank2=observational_identifiability(np.array([[1.,1.],[2.,2.],[3.,3.],[1.,0.]]),np.ones(2),np.ones(4))
    assert min(entropy)>=-1e-12 and max(defects)<1e-10 and max(mass)<1e-12
    assert max(gauge)<1e-10 and equilibrium_flux<1e-12 and abs(cycle_logK_sum)<1e-12
    assert min(inward)>=0
    result=dict(schema='synthetic_thermodynamic_foundation_audit_v1',created_utc=datetime.now(timezone.utc).isoformat(),
        evidence_class='synthetic_mathematical_fixture_only',physical_rates_validated=False,
        experimental_observations=0,new_quantum_calculations=0,ode_integrations=0,
        species_note='A/B/C are abstract inventory-preserving labels, not proposed chemical species',
        inputs=dict(seed=20260916,temperature_K=383.15,standard_concentration_M=1.,
                    species_mu0_kcal_mol=[0.,1.,-1.],transition_mu0_kcal_mol=[22.,23.,22.5],
                    alpha=alpha.tolist(),beta=beta.tolist(),random_positive_states=200,
                    concentration_sampling_M='log-uniform 1e-6 to 1, mathematical stress test only'),
        checks=dict(max_dissipation_identity_defect_kcal_L_s=max(defects),
                    min_entropy_production_kcal_L_K_s=min(entropy),
                    max_inventory_derivative_defect_M_s=max(mass),
                    max_energy_reference_shift_rhs_defect_M_s=max(gauge),
                    max_equilibrium_net_flux_M_s=equilibrium_flux,cycle_logK_sum=cycle_logK_sum,
                    all_three_zero_concentration_boundaries_point_inward=min(inward)>=0),
        observation_examples={'indistinguishable_parallel_parameters':rank1,'additional_discriminating_observation':rank2},
        limitations=['Synthetic fixtures prove no catalytic mechanism or parameter calibration',
                     'Ideal, closed, isothermal interior dissipation only',
                     'Local SVD rank is not global structural identifiability',
                     'No model selection, experimental optimization or real TOF performed'])
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    args.output.mkdir(parents=True)
    sources=[Path(__file__),REPO/'src/pincer_catmech/kinetics/thermodynamic_network.py']
    result['source_sha256']={p.relative_to(REPO).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,default=convert,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2,default=convert,allow_nan=False))


if __name__=='__main__':
    main()
