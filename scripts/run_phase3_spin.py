"""Sequential, evidence-preserving vertical spin campaign; never fabricates MECPs.

Example (campaign runtime with NumPy/ASE/psutil):
  python scripts/run_phase3_spin.py --python C:/path/to/phase7/python.exe
    --catalyst Fe_macho_pnp_iPr --timeout-seconds 120 --basis sto-3g

Default scope: 24 rows, 18 Fe/Co/Mn targets x 3 multiplicities; Ru records are
explicitly outside the spin-screen scope. Missing or invalid geometries remain
missing, and a failed state prevents its gap from becoming an accepted label.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

import numpy as np
import psutil
from ase.io import read
from ase.units import Bohr

from pincer_catmech.quantum.spin_provider import validate_spin_system
from pincer_catmech.quantum.identity import catalyst_identity_screen


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


@lru_cache(maxsize=4)
def engine_fingerprint(prefix):
    prefix=Path(prefix)
    paths=list((prefix/"conda-meta").glob("psi4-*.json"))
    paths+=list((prefix/"Lib/site-packages/psi4").glob("core*.pyd"))
    return {str(p.relative_to(prefix)):digest(p) for p in sorted(paths)}


def accepted_cache(path, signature):
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("signature_sha256") != signature or data.get("status") != "converged":
        return None
    for item in data.get("native_evidence", []):
        if not Path(item["path"]).is_file() or digest(item["path"]) != item["sha256"]:
            return None
    if not data.get("native_evidence"):
        return None
    responses=[item for item in data["native_evidence"] if Path(item["path"]).name=="response.json"]
    if len(responses)!=1:
        return None
    native=json.loads(Path(responses[0]["path"]).read_text(encoding="utf-8"))
    state=native["states"][0]
    if (state["energy_hartree"]!=data.get("energy_hartree")
            or state["gradient_hartree_bohr"]!=data.get("gradient_hartree_bohr")
            or state["spin_squared"]!=data.get("spin_squared")
            or native["geometry_sha256"]!=data.get("geometry_sha256")):
        return None
    return data


def state_attempt(args, catalyst_id, metadata, atoms, source, multiplicity, identity):
    symbols = atoms.get_chemical_symbols()
    positions = atoms.positions / Bohr
    protocol = dict(provider="psi4", method=args.method, basis=args.basis, solvent=None,
                    scf_type="df", e_convergence=args.e_convergence, d_convergence=args.d_convergence,
                    scf_max_iterations=args.scf_max_iterations, grid_radial=args.grid_radial,
                    grid_spherical=args.grid_spherical, scf_algorithm=args.scf_algorithm,
                    soscf_start_convergence=args.soscf_start_convergence if args.scf_algorithm=="soscf" else None,
                    soscf_max_iter=5 if args.scf_algorithm=="soscf" else None,
                    engine_fingerprint=engine_fingerprint(str(args.python.parent)),
                    reference=("rks" if multiplicity == 1 else "uks") if args.method.lower() not in ("hf", "scf")
                    else ("rhf" if multiplicity == 1 else "uhf"),
                    executable_sha256=digest(args.python), provider_source_sha256=digest(
                        Path(__file__).resolve().parents[1]/"src/pincer_catmech/quantum/spin_provider.py"))
    paired_protocol={k:v for k,v in protocol.items() if k!="reference"}
    paired_protocol_sha=hashlib.sha256(json.dumps(paired_protocol,sort_keys=True).encode()).hexdigest()
    geometry_sha = hashlib.sha256(np.asarray(positions, dtype="<f8").tobytes()+" ".join(symbols).encode()).hexdigest()
    signature = hashlib.sha256(json.dumps(dict(symbols=symbols, positions_bohr=positions.tolist(),
                    charge=metadata["charge"], multiplicity=multiplicity, protocol=protocol), sort_keys=True).encode()).hexdigest()
    canonical = args.output_dir / catalyst_id / f"multiplicity_{multiplicity}" / signature
    cached = accepted_cache(canonical / "result.json", signature)
    if cached:
        cached["reused_verified_cache"] = True
        return cached
    attempt = canonical / f"attempt_{uuid.uuid4().hex}"
    attempt.mkdir(parents=True, exist_ok=False)
    scratch = args.scratch_dir / attempt.name
    scratch.mkdir(parents=True, exist_ok=False)
    request = dict(symbols=symbols, positions_bohr=positions.tolist(), charge=metadata["charge"],
                   multiplicities=[multiplicity], single_state=True, method=args.method, basis=args.basis,
                   memory_mib=args.memory_mib, threads=2, scf_max_iterations=args.scf_max_iterations,
                   grid_radial=args.grid_radial, grid_spherical=args.grid_spherical,
                   scf_algorithm=args.scf_algorithm,
                   soscf_start_convergence=args.soscf_start_convergence,
                   e_convergence=args.e_convergence,d_convergence=args.d_convergence,
                   workdir=str(attempt / "native"))
    dump(scratch / "request.json", request)
    command = [str(args.python), "-m", "pincer_catmech.quantum.spin_provider", "--request",
               str(scratch/"request.json"), "--result", str(attempt/"response.json")]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
               PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    prefix = args.python.parent
    env["PATH"] = str(prefix/"Library/bin") + os.pathsep + str(prefix) + os.pathsep + env.get("PATH", "")
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]/"src")
    status, failure, peak_rss = "failed", None, 0
    start = time.monotonic()
    started_utc=datetime.now(timezone.utc).isoformat()
    with (attempt/"worker.out").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(command, cwd=scratch, env=env, stdout=log, stderr=subprocess.STDOUT)
        measured = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                peak_rss = max(peak_rss, measured.memory_info().rss)
            except psutil.NoSuchProcess:
                break
            if time.monotonic()-start >= args.timeout_seconds:
                status, failure = "timeout", "Per-state wall-time limit"
                process.kill()
                break
            if peak_rss > args.rss_limit_mib * 1024**2:
                status, failure = "memory_limit", "Observed resident memory exceeded pilot limit"
                process.kill()
                break
            time.sleep(.1)
        code = process.wait()
    row = dict(catalyst_id=catalyst_id, family=metadata.get("backbone"), metal=metadata.get("metal"),
               chemical_state=args.chemical_state, multiplicity=multiplicity, charge=metadata["charge"],
               signature_sha256=signature, geometry_sha256=geometry_sha, source_xyz=str(source),
               source_xyz_sha256=digest(source), symbols=symbols, positions_bohr=positions.tolist(),
               source_metadata_path=str(source.parent/"metadata.json"),
               source_metadata_sha256=digest(source.parent/"metadata.json"),
               protocol=protocol, protocol_sha256=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest(),
               paired_protocol_sha256=paired_protocol_sha,
               identity_validation=identity, status=status, failure=failure, returncode=code,
               intended_catalyst_identity_pass=bool(identity.get("retained",False)),
               wall_seconds=time.monotonic()-start, peak_rss_mib=peak_rss/1024**2,
               memory_configured_mib=args.memory_mib, timeout_seconds=args.timeout_seconds,
               started_utc=started_utc, native_directory=str(attempt),
               native_evidence=[], reused_verified_cache=False,
               evidence="vertical electronic spin diagnostic; no free-energy barrier, MECP, or experiment")
    if code == 0 and (attempt/"response.json").is_file():
        response = json.loads((attempt/"response.json").read_text(encoding="utf-8"))
        state = response["states"][0]
        if response["geometry_sha256"] != geometry_sha:
            row["failure"] = "Worker geometry hash mismatch"
        else:
            row.update(state)
            row["status"] = "converged"
            row["psi4_version"] = response["psi4_version"]
            row["spin_contamination"] = state["spin_squared"]-state["expected_spin_squared"]
            row["spin_contamination_flag"] = abs(row["spin_contamination"]) > .1
            row["wall_seconds"] = time.monotonic()-start
    elif failure is None:
        failed = list((attempt/"native").glob("pair_*/failure.json"))
        if failed:
            detail = json.loads(failed[0].read_text(encoding="utf-8"))
            row["failure"] = detail["error"]
            row["status"] = "scf_not_converged" if detail["error_type"] == "SCFConvergenceError" else "backend_failure"
        else:
            row["failure"] = "Worker failed before a native result was available"
            row["status"] = "backend_failure"
    for path in sorted(attempt.rglob("*")):
        if path.is_file() and path.name in {"psi4.out", "worker.out", "response.json", "request.json", "failure.json"}:
            row["native_evidence"].append(dict(path=str(path),sha256=digest(path)))
    dump(attempt / "record.json", row)
    # Canonical index can be updated, but every attempted calculation remains immutable.
    dump(canonical / "result.json", row)
    return row


def summarize(args, rows, excluded):
    gaps=[]
    for catalyst in sorted({r["catalyst_id"] for r in rows}):
        records = [r for r in rows if r["catalyst_id"]==catalyst]
        for low, high in ((1,3),(1,5),(3,5)):
            a = next((r for r in records if r["multiplicity"]==low),None)
            b = next((r for r in records if r["multiplicity"]==high),None)
            valid = bool(a and b and a["status"]==b["status"]=="converged"
                         and a["geometry_sha256"]==b["geometry_sha256"] and a["charge"]==b["charge"]
                         and a.get("paired_protocol_sha256")
                         and a["paired_protocol_sha256"]==b.get("paired_protocol_sha256"))
            item=dict(catalyst_id=catalyst,multiplicity_low=low,multiplicity_high=high,
                      status="computed_vertical_gap" if valid else "missing_or_failed_state",
                      gap_high_minus_low_hartree=None, mecp_status="not_searched",
                      diagnostic_label_eligible=False, validated_catalyst_label=False)
            if valid:
                item.update(gap_high_minus_low_hartree=b["energy_hartree"]-a["energy_hartree"],
                            geometry_sha256=a["geometry_sha256"], charge=a["charge"],
                            source_state_signatures=[a["signature_sha256"],b["signature_sha256"]],
                            spin_contamination_flag=a["spin_contamination_flag"] or b["spin_contamination_flag"],
                            intended_catalyst_identity_pass=a["intended_catalyst_identity_pass"] and b["intended_catalyst_identity_pass"],
                            diagnostic_label_eligible=(not a["spin_contamination_flag"] and not b["spin_contamination_flag"]
                                and a["intended_catalyst_identity_pass"] and b["intended_catalyst_identity_pass"]),
                            paired_protocol_sha256=a["paired_protocol_sha256"],
                            evidence="DFT vertical electronic gap at the source geometry, not a MECP gap or barrier")
            gaps.append(item)
    summary=dict(method=args.method,basis=args.basis,chemical_state=args.chemical_state,
                 grid_radial=args.grid_radial,grid_spherical=args.grid_spherical,scf_algorithm=args.scf_algorithm,
                 soscf_start_convergence=args.soscf_start_convergence if args.scf_algorithm=="soscf" else None,
                 e_convergence=args.e_convergence,d_convergence=args.d_convergence,
                 target_states=len(rows),states_converged=sum(r["status"]=="converged" for r in rows),
                 planned_target_states=args.planned_states,run_complete=len(rows)==args.planned_states,
                 states_failed_or_missing=sum(r["status"]!="converged" for r in rows),
                 state_records=rows,spin_gaps=gaps,excluded=excluded,
                 mecp_converged=0,quantitative_spin_ranking_validated=False,
                 interpretation="Bounded small-basis vertical diagnostic; failed states cannot establish absent crossings")
    dump(args.output_dir/"summary.json",summary)
    fields=["catalyst_id","multiplicity","charge","status","energy_hartree","spin_squared","spin_contamination",
            "spin_contamination_flag","wall_seconds","peak_rss_mib","source_xyz_sha256","signature_sha256","failure"]
    with (args.output_dir/"spin_states.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    with (args.output_dir/"vertical_gaps.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["catalyst_id","multiplicity_low","multiplicity_high","status",
                             "gap_high_minus_low_hartree","spin_contamination_flag","mecp_status"],extrasaction="ignore")
        writer.writeheader(); writer.writerows(gaps)
    if args.aggregate_path.resolve() != (args.output_dir/"summary.json").resolve():
        aggregate_results(args.aggregate_path,args.output_dir/"summary.json")
    return summary


def aggregate_results(path, active_matrix_path):
    """Keep one stable entry point without flattening distinct QC protocols."""
    base=Path(path).parent
    matrix=json.loads(Path(active_matrix_path).read_text(encoding="utf-8"))
    rows=matrix["state_records"]
    diagnostics=[]
    for source in sorted(base.glob("*/result.json")):
        record=json.loads(source.read_text(encoding="utf-8"))
        diagnostics.append(dict(path=str(source),sha256=digest(source),
            system=record.get("system",record.get("scope")),method=record.get("method"),
            status=record.get("status","capability_diagnostic"),converged=record.get("converged",False),
            evaluations=record.get("evaluations"),minimum_verified=record.get("minimum_verified",False)))
    other_protocols=[]
    for source in sorted(base.glob("*/summary.json")):
        record=json.loads(source.read_text(encoding="utf-8"))
        other_protocols.append(dict(path=str(source),method=record.get("method"),basis=record.get("basis"),
            grid_radial=record.get("grid_radial"),grid_spherical=record.get("grid_spherical"),
            scf_algorithm=record.get("scf_algorithm"),states_converged=record.get("states_converged"),
            target_states=record.get("target_states"),sha256=digest(source)))
    data=dict(active_matrix_summary=str(active_matrix_path),active_matrix_summary_sha256=digest(active_matrix_path),
              run_complete=matrix["run_complete"],planned_target_states=matrix["planned_target_states"],
              pending_target_states=max(0,matrix["planned_target_states"]-len(rows)),
              target_states_recorded=len(rows),physical_state_attempts=sum("returncode" in r for r in rows),
              states_converged=matrix["states_converged"],states_failed_or_missing=matrix["states_failed_or_missing"],
              missing_geometry_slots=sum(r["status"]=="missing_source_geometry" for r in rows),
              s2_and_identity_eligible_states=sum(r["status"]=="converged" and not r.get("spin_contamination_flag",True)
                                                and r.get("intended_catalyst_identity_pass",False) for r in rows),
              computed_vertical_gaps=sum(g["status"]=="computed_vertical_gap" for g in matrix["spin_gaps"]),
              eligible_diagnostic_gaps=sum(g["diagnostic_label_eligible"] for g in matrix["spin_gaps"]),
              mecp_accepted=0,validated_catalyst_gaps=0,diagnostics=diagnostics,target_protocol_runs=other_protocols,
              interpretation="Distinct protocols retained separately. Engine examples are not pincer labels; vertical gaps are not MECPs.")
    dump(path,data)
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--provider",choices=["psi4"],default="psi4")
    parser.add_argument("--python","--psi4-python",dest="python",type=Path,required=True,help="Existing Psi4 Python interpreter")
    parser.add_argument("--structures-dir",type=Path,default=Path("data/structures"))
    parser.add_argument("--output-dir",type=Path,default=Path("data/phase3/spin/pincer_vertical_pbe_sto3g_001"))
    parser.add_argument("--scratch-dir",type=Path,default=Path("work/phase3/spin/target_scratch"))
    parser.add_argument("--aggregate-path",type=Path,default=Path("data/phase3/spin/summary.json"))
    parser.add_argument("--method",default="pbe")
    parser.add_argument("--basis",default="sto-3g")
    parser.add_argument("--chemical-state",default="active",choices=["active","hydrogenated","protonated_reference"])
    parser.add_argument("--timeout-seconds",type=float,default=120)
    parser.add_argument("--memory-mib",type=int,default=500)
    parser.add_argument("--rss-limit-mib",type=float,default=600)
    parser.add_argument("--scf-max-iterations",type=int,default=120)
    parser.add_argument("--grid-radial",type=int,default=75)
    parser.add_argument("--grid-spherical",type=int,default=194)
    parser.add_argument("--scf-algorithm",choices=["diis","soscf"],default="diis")
    parser.add_argument("--soscf-start-convergence",type=float,default=1e-4)
    parser.add_argument("--e-convergence",type=float,default=1e-9)
    parser.add_argument("--d-convergence",type=float,default=1e-8)
    parser.add_argument("--catalyst",action="append",help="Specific catalyst(s) for measured pilot before full matrix")
    parser.add_argument("--first-catalyst",help="Evaluate this measured pilot target first while retaining full selected matrix")
    parser.add_argument("--multiplicities",type=int,nargs="+",default=[1,3,5])
    args=parser.parse_args()
    for name in ("python","structures_dir","output_dir","scratch_dir","aggregate_path"):
        setattr(args,name,getattr(args,name).resolve())
    if not args.python.is_file() or not args.structures_dir.is_dir():
        parser.error("Existing Python and structure directory required")
    if args.timeout_seconds<=0 or not 250<=args.memory_mib<=600 or not 100<=args.rss_limit_mib<=600:
        parser.error("Positive timeout, 250-600 MiB configuration and 100-600 MiB RSS limit required")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    rows,excluded=[],[]
    candidates=sorted(p for p in args.structures_dir.iterdir() if p.is_dir()
                      and (args.catalyst is None or p.name in args.catalyst))
    if args.catalyst and set(args.catalyst)-{p.name for p in candidates}:
        parser.error("Requested catalyst identifiers were not found")
    if args.first_catalyst:
        if args.first_catalyst not in {p.name for p in candidates}:
            parser.error("First catalyst is outside the selected matrix")
        candidates.sort(key=lambda p:(p.name!=args.first_catalyst,p.name))
    args.planned_states=sum(p.name.startswith(("Fe_","Co_","Mn_")) for p in candidates)*len(args.multiplicities)
    for candidate in candidates:
        if not candidate.name.startswith(("Fe_","Co_","Mn_")):
            excluded.append(dict(catalyst_id=candidate.name,status="outside_Fe_Co_Mn_spin_scope")); continue
        folder=candidate/args.chemical_state
        source,meta_path=folder/"best_found.xyz",folder/"metadata.json"
        if not source.is_file() or not meta_path.is_file():
            for multiplicity in args.multiplicities:
                rows.append(dict(catalyst_id=candidate.name,multiplicity=multiplicity,status="missing_source_geometry",
                                 failure=f"Missing {args.chemical_state} best_found.xyz or metadata.json"))
            summarize(args,rows,excluded)
            continue
        metadata=json.loads(meta_path.read_text(encoding="utf-8"))
        atoms=read(source)
        identity=catalyst_identity_screen(atoms,metadata)
        # Identity diagnostics are retained. A false coordination/bond graph
        # is not silently promoted to the intended pincer structure.
        for multiplicity in args.multiplicities:
            try:
                validate_spin_system(atoms.get_chemical_symbols(),atoms.positions/Bohr,
                                     metadata["charge"],[multiplicity],require_pair=False)
                row=state_attempt(args,candidate.name,metadata,atoms,source,multiplicity,identity)
            except (ValueError,OSError,RuntimeError) as exc:
                row=dict(catalyst_id=candidate.name,multiplicity=multiplicity,status="invalid_input_or_launch_failure",failure=str(exc))
            rows.append(row)
            summarize(args,rows,excluded)
            print(json.dumps({k:row.get(k) for k in ("catalyst_id","multiplicity","status","energy_hartree",
                             "wall_seconds","peak_rss_mib","failure")}),flush=True)
    summary=summarize(args,rows,excluded)
    print(json.dumps({k:summary[k] for k in ("target_states","states_converged","states_failed_or_missing")}),flush=True)
    return 0


if __name__=="__main__":
    sys.exit(main())
