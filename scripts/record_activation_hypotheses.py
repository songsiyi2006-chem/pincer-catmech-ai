"""Balanced activation hypotheses for the public precursor, not computed paths."""
from pathlib import Path
from collections import Counter
import hashlib
import json

REPO=Path(__file__).resolve().parents[1]


def main():
    target=REPO/'data/phase4/activation_hypotheses.json'
    if target.exists():
        raise FileExistsError('Hypothesis ledger is immutable; choose a new version')
    precursor={'C':18,'H':37,'Br':1,'Mn':1,'N':1,'O':2,'P':2}
    deprotonated=dict(precursor,H=36)
    br_lost={k:v for k,v in deprotonated.items() if k!='Br'}
    species={
        'P0':{'composition':precursor,'charge':0,'status':'public_crystal_precursor_charge_working_assignment'},
        'KOtBu_formula_unit':{'composition':{'K':1,'C':4,'H':9,'O':1},'charge':0,'status':'stoichiometric_reagent_unit_not_solution_aggregation_claim'},
        'K_plus':{'composition':{'K':1},'charge':1,'status':'hypothetical_explicit_counterion'},
        'P_minus_H_with_Br':{'composition':deprotonated,'charge':-1,'status':'hypothetical_deprotonated_Br_retained_state'},
        'A0_Br_lost':{'composition':br_lost,'charge':0,'status':'hypothetical_deprotonated_Br_dissociated_state'},
        'KBr_formula_unit':{'composition':{'K':1,'Br':1},'charge':0,'status':'stoichiometric_salt_ledger_phase_and_activity_undetermined'},
        'tBuOH':{'composition':{'C':4,'H':10,'O':1},'charge':0,'status':'neutral_alcohol_ledger'},
    }
    reactions=[
        dict(id='H1_Br_retention',reactants={'P0':1,'KOtBu_formula_unit':1},
             products={'K_plus':1,'P_minus_H_with_Br':1,'tBuOH':1},
             note='Explicit ions specify stoichiometry only; free, paired and aggregated forms require distinct actual models'),
        dict(id='H2_Br_dissociation',reactants={'P0':1,'KOtBu_formula_unit':1},
             products={'A0_Br_lost':1,'KBr_formula_unit':1,'tBuOH':1},
             note='Br release and KBr phase must be verified; neutral activated complex is not obtained by deleting H alone')]
    def totals(side):
        atoms=Counter();charge=0
        for name,n in side.items():
            atoms.update({e:n*v for e,v in species[name]['composition'].items()})
            charge+=n*species[name]['charge']
        return dict(atoms),charge
    for r in reactions:
        left,right=totals(r['reactants']),totals(r['products'])
        assert left==right
        r.update(element_and_charge_balanced=True,balanced_atoms=left[0],balanced_charge=left[1],
                 thermodynamic_or_kinetic_parameters=None,validated_activation_mechanism=False)
    source=REPO/'data/phase4/public_structure/import_v002/Mn1/metadata.json'
    source_metadata=json.loads(source.read_text(encoding='utf-8'))
    if (source_metadata.get('formula')!='C18H37BrMnNO2P2' or source_metadata.get('charge')!=0
            or source_metadata.get('proton_site_element')!='N'):
        raise ValueError('Source is not the declared neutral NH-bearing precursor')
    data=dict(schema='precursor_activation_hypothesis_ledger_v1',evidence_class='balanced_hypotheses_only',
              source_metadata=source.relative_to(REPO).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              species=species,reactions=reactions,new_quantum_jobs=0,accepted_activation_paths=0,
              limitations=['Integer balance does not prove existence, kinetics, oxidation state or solution speciation',
                  'No isolated proton reservoir or aqueous pH is assumed',
                  'Ion pairing, aggregation, solubility and KBr solid activity remain unresolved',
                  'Both hypotheses can be wrong; qualified chemistry review precedes path generation'],
              implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'balanced_hypotheses':2,'accepted_activation_paths':0,'new_quantum_jobs':0}))


if __name__=='__main__':
    main()
