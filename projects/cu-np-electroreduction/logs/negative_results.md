# Retained failures and limitations / 保留的失败与限制

- **EXECUTED RESULT:** First `run_pilot.py` preparation failed with
  `ValueError: Invalid suffix 'provenance.json'` before a native calculation.
  Fixed to `.provenance.json`; deterministic structure generation reran. No
  energy or successful quantum job is attributed to that failed preparation.
- **EXECUTED RESULT:** An initial unittest process exited 1 while entering the
  NumPy-backed three-state steady-state check, without a Python exception.
  Root cause was not established. Replaced this unnecessary dependency with
  the analytical three-state spanning-tree expression; independent balance
  and equilibrium tests validate the expression. No environment was modified.
- **EXECUTED RESULT:** P000 optimized but P009 had −50.49 cm−1: minimum rejected.
  The initial parser omitted the six unlabelled rigid-body modes and compared
  54 modes with an incorrect count of 60. It also rejected the geometry; the
  real imaginary frequency independently requires rejection. The corrected
  parser and regenerated summary retain all 60 modes, and raw files remain
  unchanged. This was not a transition-state search or a validated TS.
- **EXECUTED RESULT:** R000 +60° methyl seed and verytight optimization, followed
  by R001 Hessian, pass the stated minimum and connectivity checks. This is a
  minimum on one gas-phase GFN2 surface, not experimental prevalence.
- **EXECUTED RESULT:** GFN1/GFN2 charging-response spread exceeds the 0.20 eV
  screening gate. Do not use these values as redox potentials or catalyst ranks.
- **EXECUTED RESULT:** Cu original full text/SI remain inaccessible; the author
  scheme supplies partial conditions only. No published catalyst reproduction,
  Cu ensemble, TS/IRC, physical microkinetics, ML or institutional HPC run occurred.
- **MODEL INFERENCE:** The full research success definition is not yet met.
  These are evidence-audit and local methodological preflight deliverables.
