"""Regenerate evidence summary from native outputs, preserving failed checks."""
import json
import math
import re
from pathlib import Path
from science import HARTREE_EV,parse_xtb,sha256
from run_pilot import ROOT,dump

def main():
    cfg=json.loads((ROOT/"config/pilot.json").read_text())
    rows=[];native_index=[]
    for record_path in sorted((ROOT/"data/raw").glob("pilot_*/P*/record.json"))+sorted((ROOT/"data/raw").glob("pilot_*/R*/record.json")):
        record=json.loads(record_path.read_text(encoding="utf-8"));folder=record_path.parent
        for name,digest in record["output_hashes"].items():
            if sha256(folder/name)!=digest:raise RuntimeError(f"Output changed: {folder/name}")
        missing_logs=[name for name in ("stdout.log","stderr.log") if not (folder/name).exists()]
        raw="\n".join((folder/name).read_text(encoding="utf-8",errors="replace") if (folder/name).exists() else "" for name in ("stdout.log","stderr.log"))
        row={"id":record["calculation_id"],"batch":folder.parent.name,"settings":record["settings"],"energy_hartree":None,"elapsed_seconds":record["elapsed_seconds"],"input_sha256":record["identity"]["input_sha256"],"status":record["status"]}
        if missing_logs: row["file_missing"]=missing_logs
        try:
            parsed=parse_xtb(raw,record["settings"]["task"])
            if record["status"]=="completed" and not missing_logs: row["energy_hartree"]=parsed["energy_hartree"]
            elif record["status"]=="completed": row["status"]="parse_rejected"
        except ValueError as exc:
            row["parse_error"]=str(exc)
            if record["status"]=="completed": row["status"]="parse_rejected"
        if record["settings"]["task"]=="hess" and (folder/"vibspectrum").exists():
            freq=[float(x) for x in re.findall(r"^\s*\d+\s+(?:[ai]\s+)?(-?\d+\.\d+)", (folder/"vibspectrum").read_text(encoding="utf-8",errors="replace"),re.M)]
            vib=[float(x) for x in re.findall(r"^\s*\d+\s+[ai]\s+(-?\d+\.\d+)", (folder/"vibspectrum").read_text(encoding="utf-8",errors="replace"),re.M)]
            min_threshold=cfg["predeclared_tolerances"]["neutral_minimum_imaginary_cm1_threshold"]
            row.update(frequency_count=len(freq),minimum_frequency_cm1=min(freq) if freq else None,vibrational_mode_count=len(vib),minimum_vibrational_frequency_cm1=min(vib) if vib else None,significant_imaginary_modes=sum(x < min_threshold for x in freq),minimum_accepted=(len(freq)==60 and len(vib)==54 and min(freq)>=min_threshold))
        rows.append(row)
        native_index.append({"id":row["id"],"record":str(record_path.relative_to(ROOT)).replace("\\","/"),"record_sha256":sha256(record_path),"files":record["output_hashes"]})
    valid_minimum=any(r["batch"]=="pilot_002" and r["status"]=="completed" and r["energy_hartree"] is not None and r.get("minimum_accepted",False) for r in rows)
    revised=[r for r in rows if valid_minimum and r["batch"]=="pilot_002" and r["settings"]["task"]=="sp" and r["energy_hartree"] is not None and r["status"]=="completed"]
    charging=[]
    for method in (1,2):
        for solvent in ("gas","dmf"):
            pair={r["settings"]["charge"]:r for r in revised if r["settings"]["gfn"]==method and r["settings"]["solvent"]==solvent}
            if len(pair)!=2: continue
            assert pair[0]["input_sha256"]==pair[-1]["input_sha256"]
            charging.append({"gfn":method,"solvent":solvent,"delta_E_qminus1_minus_q0_eV":(pair[-1]["energy_hartree"]-pair[0]["energy_hartree"])*HARTREE_EV,"calculation_ids":[pair[-1]["id"],pair[0]["id"]],"label":"EXECUTED RESULT"})
    spread={s:abs(next(x for x in charging if x["gfn"]==1 and x["solvent"]==s)["delta_E_qminus1_minus_q0_eV"]-next(x for x in charging if x["gfn"]==2 and x["solvent"]==s)["delta_E_qminus1_minus_q0_eV"]) for s in ("gas","dmf") if len([x for x in charging if x["solvent"]==s])==2}
    # Distance-based connectivity audit of all atoms against generated molecular graph.
    from rdkit import Chem
    mol=Chem.AddHs(Chem.MolFromSmiles(cfg["scaffold"]["smiles"]))
    lines=(ROOT/"data/raw/pilot_002/R000-rotor-reopt/xtbopt.xyz").read_text().splitlines()[2:]
    coords=[[float(v) for v in x.split()[1:4]] for x in lines]
    expected={tuple(sorted([b.GetBeginAtomIdx(),b.GetEndAtomIdx()])) for b in mol.GetBonds()}
    table=Chem.GetPeriodicTable();found=set()
    for i in range(mol.GetNumAtoms()):
        for j in range(i):
            distance=math.sqrt(sum((coords[i][k]-coords[j][k])**2 for k in range(3)))
            cutoff=1.25*(table.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum())+table.GetRcovalent(mol.GetAtomWithIdx(j).GetAtomicNum()))
            if distance < cutoff: found.add((j,i))
    threshold=cfg["predeclared_tolerances"]["charging_energy_method_spread_flag_eV"]
    initial_rejected=any(r["batch"]=="pilot_001" and r.get("minimum_accepted") is False for r in rows)
    summary={"label":"EXECUTED RESULT","scope":"molecular fixed-charge preflight; not Cu-site reproduction, no kinetics or electrochemical potential", "native_calls":len(rows),"completed_native_calls":sum(r["status"]=="completed" for r in rows),"elapsed_sum_seconds":sum(r["elapsed_seconds"] for r in rows),"initial_minimum_rejected":initial_rejected,"identity_connectivity_preserved":found==expected,"connectivity_note":"distance graph only, not bond-order, chemical or active-site validation", "rows":rows,"fixed_geometry_charging_energies":charging,"hamiltonian_spread_eV":spread,"predeclared_spread_threshold_eV":threshold,"method_sensitivity_stop":any(v>threshold for v in spread.values()),"interpretation":"These are E(anion)-E(neutral) with the same nuclei. ALPB uses equilibrium implicit solvent polarization. This is NOT a nonequilibrium vertical electron affinity, electrode potential, ET barrier, free energy, selectivity, or Cu mechanism. No electron chemical potential/reference calibration is supplied."}
    dump(ROOT/"results/pilot_summary.json",summary);dump(ROOT/"data/raw/index.json",native_index)
    print(json.dumps({k:summary[k] for k in ["native_calls","elapsed_sum_seconds","identity_connectivity_preserved","fixed_geometry_charging_energies","hamiltonian_spread_eV","method_sensitivity_stop"]},indent=2))
if __name__=="__main__":main()
