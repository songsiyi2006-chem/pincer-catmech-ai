# PINCER-CATMECH-AI

Reproducible qRRHO thermochemistry and three-dimensional pincer descriptors for
Gaussian 16 / ORCA workflows. This initialization provides numerical building
blocks, not a trained model or an experimentally validated catalyst predictor.

## Install and verify

Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest tests/ -v
python examples/mock_pipeline.py
```

On Windows use `.venv\Scripts\activate` in cmd, or
`.\.venv\Scripts\Activate.ps1` in PowerShell. Git Bash can run the deployment
script and automatically discovers either venv layout.

## Thermochemistry

```python
from pincer_catmech.kinetics.free_energy import parse_thermochemistry

result = parse_thermochemistry("calculation.log")
print(result.G_298_qRRHO_sol)
```

Energy fields are kcal/mol; entropy fields are cal/(mol K). The historical
`298` names describe the default temperature; inspect `temperature` for a
calculation performed at another temperature. Read the module docstrings for
supported parser layouts, job selection and explicit transition-state policy.
Inputs need normally terminated frequency jobs and complete thermochemistry.
The parser reconstructs harmonic entropy before damping ORCA outputs, whose
reported entropy may already use qRRHO. Ambiguous/incomplete jobs raise errors.
External format checks passed on one Gaussian 16 C.01 file (57 modes) and one
ORCA 4.0.1 file (one positive and five projected external modes). ORCA 6.1
handling is based on its official documented layout and synthetic regression
tests; this is not an all-version compatibility guarantee. See
`docs/EXTERNAL_FORMAT_CHECKS.json` for pinned source URLs and hashes.

The default uses Grimme's reduced effective moment of inertia. An explicit
uncapped option reproduces the inertia expression supplied in the mission.
Changing standard state to 1 M **does not compute a solvation free energy**.
For a reaction, apply molecular corrections with stoichiometric coefficients.

## Steric features

```python
from ase.io import read
from pincer_catmech.features.steric import bite_angle, buried_volume

atoms = read("complex.xyz")
angle = bite_angle(atoms, 0, 1, 2)
profile = buried_volume(atoms, 0, ligand_indices=[1, 2, 3])
print(angle, profile.percent_buried_volume, profile.confidence_interval_percent)
```

Indices are zero-based and coordinates are angstroms. The example atom indices
must be replaced with the actual metal, donor and complete ligand selection for
your structure. The central metal is always excluded from occupied volume.
Hydrogens are omitted by default. Bondi radii are scaled by 1.17 in a 3.5 Å
sphere. Unsupported elements require explicit radii; ligand atoms and spectator
species cannot be inferred from an XYZ coordinate file alone.

## Reproducible demonstration

`examples/mock_pipeline.py` and `tests/test_pipeline.py` use deliberately synthetic
Ru-PNP-like coordinates and frequencies, including 25.4, 68.2 and 120.5 cm⁻¹.
They verify numerical plumbing without claiming a real optimized Ru complex.
Analytical tests independently check entropy limits, standard-state shifts,
geometric angles and sphere volumes. Monte Carlo intervals describe sampling
uncertainty only, not uncertainty in a molecular geometry or a radius model.

## Deployment

```bash
bash scripts/deploy.sh
```

The script requires branch `main`, synchronizes an existing `origin/main` with
fast-forward only, creates a local venv if absent, installs this package and its
test dependency, and runs `python -m pytest tests/ -v`. Only after success does
it stage the project paths, commit with the specified conventional message,
and execute `git push origin main`. It never force-pushes. An upstream change
that conflicts with local work or a failed test aborts the deployment.
Use `PYTHON` to select a bootstrap interpreter and `PINCER_VENV` to select a
different venv directory. Existing environments are reused, not recreated.

See [English report](docs/TECHNICAL_REPORT_EN.md) and
[中文技术报告](docs/TECHNICAL_REPORT_ZH.md) for derivations, citations, algorithms,
evidence limits and the intended connection to pincer catalysis research.
