# Cu–N–P organic electroreduction / Cu–N–P 有机电还原

**2026-09-20: evidence audit and local molecular preflight, not a completed
catalytic mechanism. 文献审计与本地分子预检已交付，催化机理和性能预测尚未完成。**

Read the [中文技术报告](reports/TECHNICAL_REPORT_ZH.md) or
[English technical report](reports/TECHNICAL_REPORT_EN.md). The original
pincer project remains separate; these are additive research files.

## Current evidence

- **LITERATURE FACT:** The target is quinazolinone–aryl-nitrile reductive
  coupling/ring opening, not an isolated two-electron hydrogenation. Partial
  source conditions are DMF/TEOA/TBAPF6, 10 mA, carbon anode/CuN3P1 cathode.
- **EXECUTED RESULT:** Five required DOI anchors plus later original studies
  are audited in a [12-study matrix](literature/evidence_matrix.md). Only four
  studies meet the stricter organic-electroreduction category (shortfall four
  to eight); eight count only with mechanism/method neighbors. Just one is the
  same quinazolinone system. Cu full text/SI/model coordinates remain unavailable.
- **EXECUTED RESULT:** Twenty native xTB calls on independently generated
  N3-methylquinazolinone: ten initial calls, rejected −50.49 cm−1 minimum,
  then ten repair/recheck calls. The repaired neutral reference has 54 positive
  vibrations (lowest 78.64 cm−1). [Raw index](data/raw/index.json).
- **EXECUTED RESULT:** GFN1/GFN2 fixed-nuclei charging differences disagree by
  0.78–0.83 eV. This exceeds the 0.20 eV gate; these values are not calibrated
  reduction potentials, Cu-site effects or validated performance predictions.
- **EXECUTED RESULT:** 31 meaningful scientific/software tests pass.
  [Independent review](reports/CODE_REVIEW.md), [test log](logs/tests.txt).
- **PLANNED CALCULATION:** Published-model reproduction, Cu ensembles,
  interface kinetics, prospective prediction and ML remain gated. Server
  access, allocation, software licenses and experimental cooperation are unverified.

## Deliverables in contractual order

| # | Deliverable | Location and actual status |
|---|---|---|
| 1 | Project charter | Report §1: narrow question, competing hypotheses, application and gates |
| 2 | Evidence matrix | `literature/`: verified sources, access gaps, role/source boundaries |
| 3 | Model specification | Report §3; `config/pilot.json`, `config/production_gate.json` |
| 4 | Pilot execution | Report §4; native `data/raw/pilot_001` and `pilot_002`, portable scripts |
| 5 | Mechanism comparison | Report §5 and `reports/DESIGN_PROCESS_*`: competing models, none selected |
| 6 | Predictive results | Report §6: quantitative predictions withheld; `config/prospective_holdout.json` frozen |
| 7 | Experimental discrimination | Report §7: five high-value experiments, controls and contrary outcomes |
| 8 | Industrial assessment | Report §8; 18 explicit scenarios, unknown electron count/product mass retained |
| 9 | Reproducibility | Report §9; environment, native input/output records, hashes, commands and tests |
| 10 | Limitations/negative results | Report §10; `logs/negative_results.md`, failed scientific acceptance preserved |

## Regenerate without running new quantum chemistry

Use the existing RDKit environment. On the audited Windows host:

```powershell
$py = 'C:\Users\HUIWEI\.codex\tools\chem-ai4s\venv\Scripts\python.exe'
# Set location to this projects/cu-np-electroreduction directory first.
& $py scripts/summarize_pilot.py
& $py scripts/process_scenarios.py
& $py -m unittest discover -s tests -v
& $py scripts/verify_package.py
```

On a configured Linux machine use `python` in place of `& $py`. The archived
summary intentionally summarizes the frozen `pilot_001` and `pilot_002`
experiments. No synthetic fixtures are loaded as physical data. The kinetic
diagnostic interface rejects physical use; read [its scope](reports/KINETICS_SCOPE.md).

## Run a fresh bounded native pilot

Do not overwrite archived records. xTB must already be installed/authorized;
Python needs RDKit (recorded version in `environment/requirements-pilot.txt`).

```powershell
$xtb = 'C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin\xtb.exe'
& $py scripts/run_pilot.py --xtb $xtb --output data/raw/replay_001
& $py scripts/refine_pilot.py --xtb $xtb --source data/raw/replay_001/P000-neutral-gfn2-gas-opt/xtbopt.xyz --output data/raw/replay_002
```

The repair is specific to this atom ordering/methyl rotor, not a general
transition-state or minimum-search algorithm. Fresh runs have their own
`index.json`; the historical summary does not automatically combine them.
Calculation identity includes input, native settings, engine and source
hashes. An unchanged existing ID verifies archived bytes and skips execution;
changed code/settings require a new output directory. Historical runner/parser
source bytes matching recorded hashes are retained under
`workflows/executed_sources/`. Use current scripts for new work; archived
source copies include the documented old Hessian-count/software limitation.

`workflows/probe_server.sh` is read-only. `workflows/pilot.slurm` is a small
molecular pilot example requiring actual paths and an approved allocation.
It is **not** a Cu periodic input or evidence of server access. A source-matched
Cu input cannot be prepared without the missing model; the next minimal B01–B04
batch is explicitly recorded as unexecuted in `config/production_gate.json`.

## Structure and evidence policy

`environment/`, `config/`, `literature/`, `structures/`, `inputs/`, `workflows/`,
`scripts/`, `tests/`, `data/raw/`, `data/processed/`, `results/`, `reports/`, `logs/`
implement the research contract. Each native record has inputs, provenance,
resources, engine version, status, convergence, parser history and downstream
links. `results/FILE_MANIFEST.json` hashes the delivered files (excluding itself,
temporary caches and later reruns); raw records additionally hash native bytes.

Evidence labels: `LITERATURE FACT`, `WORKING HYPOTHESIS`, `PLANNED CALCULATION`,
`EXECUTED RESULT`, `MODEL INFERENCE`, `REQUIRES EXPERIMENT`.
Calculations and tests do not establish experiment, catalyst accessibility,
manufacturing use, novelty or publication readiness. Third-party full texts
and figures were not republished; see [literature access policy](literature/README.md).

**Next decision / 下一步：** obtain lawful Cu original/SI/model coordinates and
verify product/half-reaction mapping, then verify server allocation and run the
four source-matched baseline geometries. Do not rank catalysts before those
gates and method-sensitivity checks pass.
