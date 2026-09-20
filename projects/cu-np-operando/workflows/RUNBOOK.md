# Reproduction and continuation / 复现与继续研究

This is a local-execution delivery. The native molecular calculations, numerical
verification cases and future electrode studies are different evidence classes.

## Environment

Reuse a Python environment with NumPy, SciPy and psutil for the workflow and
tests. Psi4 jobs run in a separate existing Psi4 interpreter. Do not install a
new environment merely to reproduce an already available environment.

Audited host paths (override on another machine):

```powershell
$py = 'C:\Users\HUIWEI\miniconda3\python.exe'
$qc = 'C:\Users\HUIWEI\miniconda3\envs\phase7\python.exe'
```

Change directory to `projects/cu-np-operando` before the project-local commands.
On Linux, use the appropriate Python executable instead of the PowerShell form.
The base interpreter above passed the numerical suite directly. If using a
different Windows Conda environment, activate it first so its `Library/bin`
DLL directory is available. A direct unactivated phase2ff call failed inside
native NumPy linear algebra during review; registering that environment's DLL
directory resolved it without installing or changing packages.

## Non-quantum checks

```powershell
& $py -m unittest discover -s tests -v
& $py scripts/evidence_gate.py gate
& $py scripts/evidence_gate.py verify
```

The expected production-readiness result is `false`. A successful software test
does not unlock physical Cu-site kinetics or a catalyst ranking.
The gate command exits with code 2 for this scientifically blocked state.

The manifest verification checks frozen bytes. Regenerating numerical audits,
reports or figures may change them; retain a new output location for an
independent replay and compare numeric quantities with stated tolerances.
Do not regenerate the archived manifest merely to hide an unexpected mismatch.

## Fresh native jobs

Use the shipped bounded runner and its `--help`. Set a NEW output directory.
Both charge states must use the same copied geometry, basis and functional to
form a vertical attachment difference. The runner refuses overwrite, checks
available RAM and limits each child by sampled process-tree RSS and wall time.

The stopping check is sampled, not an OS-enforced hard allocation limit.
The Psi4 memory setting is not total RSS. Sequential execution is intentional.
The neutral structure is a previous GFN2-xTB minimum, not an optimized DFT
minimum and not a published Cu-site geometry. Preserve failed jobs too.

## Why there is no executable Cu production input

Source-matched CuN4/CuN3P1 coordinates, exact partner/product atom mapping,
non-aqueous potential calibration and a validated constant-potential engine are
missing. Producing an apparently ready periodic input would require invented
chemistry and numerical settings. No such input has been fabricated.

To continue, acquire the lawful SI and catalyst coordinates, resolve conflicting
additive acronyms, verify the atom/electron balance, configure the institution's
actual engine and scheduler, and benchmark one source-matched cell before any
screening. Hash each source and record licensed redistribution restrictions.

## Next minimal physical batch (NOT run)

1. Match and converge two source-supported bare/adsorbed Cu reference pairs.
2. Calibrate the correct solvent/reference potential and charge treatment.
3. Compute a small set of competing adsorptions and reaction paths at fixed
   electrode potentials with solvent/configuration uncertainty.
4. Infer a target-reaction kinetic network only after connectivity and electron
   balance are established; the neutral verification kernel cannot replace it.
5. Collect site density, flow geometry and experimental transport data before
   converting per-site rates to geometric current or reactor productivity.

Existing old project files remain frozen. This project does not revise their
native calculations, acceptance thresholds or reported negative results.
