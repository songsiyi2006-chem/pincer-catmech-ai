"""Explicit engineering scenarios, never experimental measurements."""
import csv
import itertools
from run_pilot import ROOT,dump
from science import electrical_kwh_kg

def main():
    rows=[]
    # Reporting per z / M is more useful than guessing unknown product chemistry.
    for voltage,fe,recovery in itertools.product([2.,4.,6.],[.3,.6,.9],[.7,.9]):
        rows.append({"label":"MODEL INFERENCE","scenario_type":"assumed_process_inputs_not_target_measurements","full_cell_V":voltage,"FE_fraction":fe,"acceptable_product_recovery":recovery,"z_electrons":"unverified","product_M_g_mol":"unverified","kWh_kg_coefficient_to_multiply_by_z_over_M":electrical_kwh_kg(voltage,1,fe,1,recovery)})
    path=ROOT/"results/process_scenarios.csv";path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    dump(ROOT/"results/process_scenario_assumptions.json",{"label":"MODEL INFERENCE","relation":"kWh/kg acceptable isolated product = coefficient * z / M(g mol^-1)","full_cell_required":True,"z":None,"M":None,"meaning":"scenario grid not measured target electricity; do not publish a target kWh/kg before balanced product identification","excludes":["cooling","pumping","solvent recovery","electrolyte recovery","separation","catalyst manufacture"],"FE_convention":"target formed; isolated recovery applied separately"})
if __name__=="__main__":main()
