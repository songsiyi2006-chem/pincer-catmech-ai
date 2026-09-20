# Cu–N–P operando research: local execution / 工况耦合研究的本地执行

**Scope agreed on 20 September 2026: complete the locally executable work and
publish it; no HPC or original SI is available. 本次完成本地可执行部分；完整恒电位
催化机制、真实排名反转、实验与工业验证尚未完成。**

Read the [中文技术报告](reports/TECHNICAL_REPORT_ZH.md) and
[English technical report](reports/TECHNICAL_REPORT_EN.md).

This is an additive continuation of the frozen
[earlier Cu molecular preflight](../cu-np-electroreduction/README.md).
Its native runs and negative results are preserved without reinterpretation.

## Evidence classes

- **CALCULATION:** new bounded Psi4 fixed-geometry gas-phase neutral/anion
  single points on the independently generated molecular substrate model.
  These are not Cu-containing structures, DFT minima or reduction potentials.
- **SOFTWARE VERIFICATION:** a reversible neutral three-state surface kernel,
  film mass-transfer coupling, conservative plug-flow integration and
  independently checked full-cell electrical accounting. Numerical fixtures do
  not represent the target reduction/coupling/ring-opening reaction.
- **LITERATURE:** updated source-access and prior-art audit, including direct
  precedents for potential-dependent ranking and multiscale coupling.
- **HYPOTHESIS:** the specific CuN4/CuN3P1 ranking boundary remains untested.
- **NOT RUN:** target periodic constant-potential paths, calibrated physical
  kinetics, trained target surrogate, real TEA/LCA, prospective experiments.

## Reproduction

See [RUNBOOK](workflows/RUNBOOK.md). Reuse installed environments; no license,
cloud allocation or external compute is assumed.

```bash
python -m unittest discover -s tests -v
python scripts/evidence_gate.py gate
python scripts/evidence_gate.py verify
```

The production gate is expected to return `physical_ranking_ready: false`.
Its CLI exit code is 2; metadata completeness cannot enable a missing target model.
This is an honest scientific stop, not a software failure. Hash validation only
binds bytes, and reviewer attestations do not turn metadata into physical proof.

See the shipped DFT runner and summary scripts for independent replay into NEW
directories. Native calculation outputs are frozen. Reconstructible integral
scratch remains on the originating workstation with an exclusion inventory;
the Git package is not a full scratch backup.

## Project map

| Location | Purpose |
|---|---|
| `reports/RESEARCH_CHARTER_ZH.md` | Original advanced requirements, not a completion certificate |
| `reports/TECHNICAL_REPORT_*.md` | Actual execution, numerical results and remaining gaps |
| `literature/` | Source evidence, novelty matrix, access log and unresolved conflicts |
| `environment/` | Live resource audit |
| `src/`, `scripts/`, `tests/` | Executable model, provenance gates and independent checks |
| `data/` | Native molecular inputs, outputs and execution records |
| `results/` | Derived audits, status, figures and byte-level manifest |
| `config/` | Unmet scientific inputs and prospective decision criteria |
| `workflows/` | Reproduction and source-dependent continuation |

Originality is not established by combining standard methods. The strongest
prior-art challenges and the narrower possible contribution are documented in
the literature audit. Software verification, electronic-structure convergence,
scientific accuracy and experimental validation are separate statements.
