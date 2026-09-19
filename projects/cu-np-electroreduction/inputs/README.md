# Inputs and execution boundary

`config/pilot.json` is the frozen initial 20-atom molecular preflight specification.
Every executed native input, command and settings file is copied into its own
`data/raw/pilot_001/P...` or `data/raw/pilot_002/R...` directory. The P batch is
retained even though its optimized geometry failed the vibrational gate.
The R batch applies a documented methyl-torsion repair and uses one identical
neutral geometry for all charge/solvent/Hamiltonian comparisons.

No published CuN3P1/CuN4 coordinate file was obtained. There is therefore no
fabricated periodic input deck here. `config/production_gate.json` lists the
missing physical inputs, minimal next batch and preregistered tolerances.
The molecular Slurm example is runnable once local paths and an approved
allocation are supplied; it is not a Cu interface production submission.
