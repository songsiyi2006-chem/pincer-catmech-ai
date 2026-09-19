"""Repair rejected methyl rotor geometry. Preserve pilot_001 without overwriting."""
import argparse
import json
import math
import os
import re
import subprocess
from pathlib import Path
from run_pilot import ROOT,execute,dump
from science import sha256

def main():
    p=argparse.ArgumentParser(); p.add_argument("--xtb",type=Path,required=True); p.add_argument("--output",type=Path,default=ROOT/"data/raw/pilot_002")
    p.add_argument("--source",type=Path,default=ROOT/"data/raw/pilot_001/P000-neutral-gfn2-gas-opt/xtbopt.xyz")
    p.add_argument("--threads",type=int,default=2);p.add_argument("--timeout",type=int,default=180);args=p.parse_args()
    args.xtb=args.xtb.resolve();args.output=args.output.resolve()
    if not 1<=args.threads<=4: p.error("Use 1-4 threads")
    source=args.source.resolve()
    lines=source.read_text().splitlines();rows=[x.split() for x in lines[2:22]]
    coords=[[float(v) for v in x[1:4]] for x in rows]
    axis=[coords[0][k]-coords[1][k] for k in range(3)]
    norm=math.sqrt(sum(x*x for x in axis));axis=[x/norm for x in axis]
    angle=math.pi/3
    # Rotate only three methyl hydrogens around N(1)-C(0), using Rodrigues.
    for idx in [12,13,14]:
        v=[coords[idx][k]-coords[0][k] for k in range(3)]
        cross=[axis[1]*v[2]-axis[2]*v[1],axis[2]*v[0]-axis[0]*v[2],axis[0]*v[1]-axis[1]*v[0]]
        dot=sum(axis[k]*v[k] for k in range(3))
        coords[idx]=[coords[0][k]+v[k]*math.cos(angle)+cross[k]*math.sin(angle)+axis[k]*dot*(1-math.cos(angle)) for k in range(3)]
    xyz=ROOT/"structures/rotor_displaced.xyz"
    content="20\nGenerated +60 degree methyl torsion after rejected -50.49 cm-1 rotor; not a literature coordinate\n"+"\n".join(rows[i][0]+" "+" ".join(f"{v:.10f}" for v in xyzrow) for i,xyzrow in enumerate(coords))+"\n"
    if xyz.exists() and xyz.read_text()!=content: raise RuntimeError("Refusing to overwrite changed seed")
    xyz.write_text(content)
    env=os.environ.copy();env["PATH"]=str(args.xtb.parent)+os.pathsep+env.get("PATH","")
    v=subprocess.run([str(args.xtb),"--version"],capture_output=True,text=True,env=env,timeout=20)
    version=re.search(r"xtb version ([^\n\r]+)",v.stdout).group(1)
    provenance={"source":str(source.relative_to(ROOT)),"source_sha256":sha256(source),"modification":"methyl hydrogen rotation +60 degrees about N-C axis", "reason":"initial Hessian -50.49 cm-1; nuclear coordinates are a local computational seed"}
    settings=dict(gfn=2,charge=0,unpaired=0,solvent="gas",task="opt",optimization="verytight")
    records=[execute("R000-rotor-reopt",xyz,settings,args,version,provenance)]
    if records[0]["status"]!="completed": raise RuntimeError("Failed repair")
    optimized=args.output/"R000-rotor-reopt/xtbopt.xyz"
    provenance={"origin_calculation_id":"R000-rotor-reopt","sha256":sha256(optimized)}
    records.append(execute("R001-rotor-hess",optimized,dict(gfn=2,charge=0,unpaired=0,solvent="gas",task="hess"),args,version,provenance))
    if records[-1].get("minimum_check")!="passed":
        dump(args.output/"index.json",{"records":records,"stop":"minimum validation failed"});raise RuntimeError("Minimum not validated; stop downstream jobs")
    n=2
    for method in (1,2):
        for solvent in ("gas","dmf"):
            for charge,unpaired in ((0,0),(-1,1)):
                settings=dict(gfn=method,charge=charge,unpaired=unpaired,solvent=solvent,task="sp")
                records.append(execute(f"R{n:03d}-gfn{method}-{solvent}-q{charge}",optimized,settings,args,version,provenance));n+=1
    dump(args.output/"index.json",{"records":records,"scope":"same-composition molecular charging sensitivity, no electrode or copper"})
    if any(x["status"]!="completed" for x in records): raise RuntimeError("At least one native job failed")
if __name__=="__main__":main()
