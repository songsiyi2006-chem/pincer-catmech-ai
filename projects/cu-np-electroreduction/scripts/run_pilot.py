"""Sequential, restartable fixed-charge molecular preflight with raw evidence."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import re
from science import parse_xtb, sha256

ROOT=Path(__file__).resolve().parents[1]
def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
def now(): return dt.datetime.now(dt.timezone.utc).isoformat()

def execute(calc_id, xyz, settings, args, version, provenance):
    folder=args.output/calc_id
    identity={"input_sha256":sha256(xyz),"settings":settings,"software_version":version,
              "runner_sha256":sha256(Path(__file__)),"parser_sha256":sha256(ROOT/"scripts/science.py")}
    if folder.exists():
        record=json.loads((folder/"record.json").read_text(encoding="utf-8"))
        if record["identity"] != identity: raise RuntimeError("Immutable calculation ID collision; choose a new output directory")
        for rel,h in record.get("output_hashes",{}).items():
            if sha256(folder/rel)!=h: raise RuntimeError("Archived output hash changed")
        return record
    folder.mkdir(parents=True)
    shutil.copyfile(xyz,folder/"input.xyz")
    command=[str(args.xtb),"input.xyz","--gfn",str(settings["gfn"]),"--chrg",str(settings["charge"]),"--uhf",str(settings["unpaired"]),"--acc","0.1","--iterations","250"]
    if settings["solvent"] != "gas": command += ["--alpb",settings["solvent"]]
    if settings["task"]=="opt": command += ["--opt",settings.get("optimization","tight")]
    elif settings["task"]=="hess": command += ["--hess"]
    record={"calculation_id":calc_id,"label":"EXECUTED RESULT","scope":"molecular_method_preflight_only","identity":identity,
        "input_structure_provenance":provenance,"command":command,"settings":settings,"resource_request":{"threads":args.threads,"timeout_seconds":args.timeout,"memory_limit":"not imposed by xTB; small 20-atom molecule"},
        "started_utc":now(),"status":"running","output_location":str(folder.relative_to(ROOT)) if folder.is_relative_to(ROOT) else str(folder),
        "downstream_results":["results/pilot_summary.json"],"parsing_history":[]}
    dump(folder/"record.json",record); dump(folder/"input_settings.json",settings)
    env=os.environ.copy(); env.update(OMP_NUM_THREADS=str(args.threads),MKL_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS=str(args.threads),OMP_STACKSIZE="256M")
    env["PATH"]=str(args.xtb.parent)+os.pathsep+env.get("PATH","")
    start=time.perf_counter()
    try:
        with (folder/"stdout.log").open("wb") as out,(folder/"stderr.log").open("wb") as err:
            proc=subprocess.run(command,cwd=folder,env=env,stdout=out,stderr=err,timeout=args.timeout)
        record["exit_code"]=proc.returncode
        text=(folder/"stdout.log").read_text(encoding="utf-8",errors="replace")+"\n"+(folder/"stderr.log").read_text(encoding="utf-8",errors="replace")
        parsed=parse_xtb(text,settings["task"])
        if proc.returncode: raise ValueError("nonzero native exit code")
        record.update(parsed); record["status"]="completed"
        if settings["task"]=="hess":
            vib=(folder/"vibspectrum").read_text(encoding="utf-8",errors="replace")
            freq=[float(x) for x in re.findall(r"^\s*\d+\s+[ai]\s+(-?\d+\.\d+)",vib,re.M)]
            record["frequencies_cm1"]=freq
            record["minimum_check"]="passed" if len(freq)==54 and min(freq)>=-20 else "not_accepted"
        record["parsing_history"].append({"utc":now(),"parser":"science.parse_xtb","result":"accepted native termination/energy"})
    except Exception as exc:
        record["status"]="failed"; record["error"]=f"{type(exc).__name__}: {exc}"
        record["parsing_history"].append({"utc":now(),"result":"failure retained"})
    record["elapsed_seconds"]=time.perf_counter()-start; record["finished_utc"]=now()
    record["output_hashes"]={p.name:sha256(p) for p in folder.iterdir() if p.is_file() and p.name!="record.json"}
    dump(folder/"record.json",record)
    print(calc_id,record["status"],round(record["elapsed_seconds"],2),flush=True)
    return record

def main():
    p=argparse.ArgumentParser(); p.add_argument("--xtb",type=Path,required=True); p.add_argument("--output",type=Path,default=ROOT/"data/raw/pilot_001")
    p.add_argument("--threads",type=int,default=2); p.add_argument("--timeout",type=int,default=180)
    args=p.parse_args(); args.xtb=args.xtb.resolve(); args.output=args.output.resolve()
    if not 1<=args.threads<=4: p.error("local preflight is bounded to 1-4 threads")
    from rdkit import Chem,rdBase
    from rdkit.Chem import AllChem,rdMolDescriptors
    cfg=json.loads((ROOT/"config/pilot.json").read_text(encoding="utf-8"))
    env=os.environ.copy(); env["PATH"]=str(args.xtb.parent)+os.pathsep+env.get("PATH","")
    v=subprocess.run([str(args.xtb),"--version"],capture_output=True,text=True,env=env,timeout=20)
    version=re.search(r"xtb version ([^\n\r]+)",v.stdout).group(1)
    xyz=ROOT/"structures/substrate_generated.xyz"
    mol=Chem.AddHs(Chem.MolFromSmiles(cfg["scaffold"]["smiles"]))
    assert rdMolDescriptors.CalcMolFormula(mol)=="C9H8N2O"
    assert sum(a.GetAtomicNum() for a in mol.GetAtoms())==84
    if not (xyz.exists() and xyz.with_suffix(".provenance.json").exists()):
        params=AllChem.ETKDGv3(); params.randomSeed=cfg["seed"]
        assert AllChem.EmbedMolecule(mol,params)==0
        assert AllChem.MMFFOptimizeMolecule(mol,maxIters=1000)==0
        xyz.parent.mkdir(parents=True,exist_ok=True); xyz.write_text(Chem.MolToXYZBlock(mol),encoding="utf-8")
        dump(xyz.with_suffix(".provenance.json"),{**cfg["scaffold"],"rdkit":rdBase.rdkitVersion,"seed":cfg["seed"],"sha256":sha256(xyz)})
    provenance=json.loads(xyz.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    if Chem.MolToSmiles(Chem.MolFromSmiles(provenance["smiles"])) != Chem.MolToSmiles(Chem.MolFromSmiles(cfg["scaffold"]["smiles"])):
        raise RuntimeError("Cached structure identity differs from requested scaffold")
    if provenance["sha256"]!=sha256(xyz): raise RuntimeError("Seed structure hash mismatch")
    opt=execute("P000-neutral-gfn2-gas-opt",xyz,{"gfn":2,"charge":0,"unpaired":0,"solvent":"gas","task":"opt"},args,version,provenance)
    records=[opt]
    if opt["status"]!="completed": raise RuntimeError("Geometry preflight failed; downstream jobs stopped")
    optimized=args.output/opt["calculation_id"]/"xtbopt.xyz"
    provenance={"origin_calculation_id":opt["calculation_id"],"sha256":sha256(optimized),"root_source":cfg["scaffold"]["source"]}
    n=1
    for method in (1,2):
        for solvent in ("gas","dmf"):
            for charge,unpaired in ((0,0),(-1,1)):
                assert (84-charge-unpaired)%2==0
                settings=dict(gfn=method,charge=charge,unpaired=unpaired,solvent=solvent,task="sp")
                records.append(execute(f"P{n:03d}-gfn{method}-{solvent}-q{charge}",optimized,settings,args,version,provenance)); n+=1
    records.append(execute("P009-neutral-gfn2-gas-hess",optimized,dict(gfn=2,charge=0,unpaired=0,solvent="gas",task="hess"),args,version,provenance))
    dump(args.output/"index.json",{"created_utc":now(),"protocol":"substrate_preflight_v1","records":records})
    if any(r["status"]!="completed" for r in records): sys.exit(2)
if __name__=="__main__": main()

