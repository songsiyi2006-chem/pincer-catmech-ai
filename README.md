# PINCER-CATMECH-AI

Reproducible qRRHO thermochemistry, pincer geometry generation, native GFN2-xTB
searches, strict saddle-point validation, and evidence-gated microkinetics.
The core initialization is followed by an actual three-hour local campaign
and a separately recorded Phase 3 extension.

**Scientific status:** the campaign produced real optimized structures and
source-verified molecular thermochemistry. The complete transition-state network,
physical TOF grid, trained barrier surrogate, and Pareto ranking remain unavailable. Software
verification is not experimental validation. The exact frozen counts and unmet
targets are recorded in [CAMPAIGN_STATUS.json](data/CAMPAIGN_STATUS.json).

Read the [English monograph](docs/MONOGRAPH_BORROWING_HYDROGEN_EN.md) or
[中文专著](docs/MONOGRAPH_BORROWING_HYDROGEN_ZH.md). The publication includes both
Markdown and typeset PDFs, actual native evidence, rejected candidates, and
portable geometry-linked datasets.

## Phase 3: spin, explicit solvent, analytical kinetics and native EGNN

Read the [English whitepaper](docs/MAGNUM_OPUS_CATMECH_EN.md),
[中文白皮书](docs/MAGNUM_OPUS_CATMECH_ZH.md), and the machine-readable
[Phase 3 status](data/phase3/PHASE3_STATUS.json). The previous campaign's frozen
status and raw evidence are preserved separately.

- The MECP implementation uses translation-invariant projected stationarity,
  an equality-constrained SQP step, damped BFGS and a merit line search. Actual
  Psi4 spin attempts and their failures are retained. Native GFN2 occupation
  changes do not supply reliable distinct-spin energy/gradient surfaces here. A small energy
  gap alone is insufficient to accept a crossing, and a crossing is not an ISC rate.
- Nine GFN2-xTB/ALPB(toluene) hemiaminal clusters optimized with one to three
  explicit tBuOH molecules; three sampled structures passed full Hessian
  minimum checks. At 383.15 K their 1 M association free energies are
  **+3.80, +10.22 and +13.69 kcal/mol**, despite attractive model association
  energies. These open hydrogen-bonded minima are not dehydration transition states.
- Two neutral hemiaminal/tBuOH dehydration seeds and one exact-band continuation
  completed 400 NEB/FIRE steps. Four original endpoint minimum certifications
  and two repeated certifications passed; all three bands failed the unchanged
  force threshold. **No transition state or activation free energy was accepted.**
  The two native archives retain all 51,322 archive-qualified member versions,
  including restart/cache files. This neutral channel does not resolve potassium
  ion pairing, free-base activity, or a cyclic proton-shuttle mechanism.
- The 18-state, 18-channel model includes a hydrogen-balanced dimer, poisoning,
  and a declared arm-cleavage hypothesis. Free base is dynamic, dimers count
  twice in the metal balance, and an eliminated cleavage-fragment ledger
  preserves atoms and charge. Forty numerical fixture grids each contain 500
  time points, exact analytic Jacobians and tangent sensitivities. Their TOF
  and control coefficients test software, not real catalyst performance.
- The native PyTorch EGNN was trained only on actual matched conformer energy
  differences: 610 structures, 540 pairs, with complete family grouping.
  Test MAE **0.19520 eV** versus **0.19400 eV** for a zero-difference baseline
  demonstrates **no predictive improvement**. Barrier and MECP-gap heads remain
  untrained. Pair predictions obey exact self-zero/antisymmetry; float64
  rotation/inversion/translation tests reach about 1.8e-15 absolute error.

Reuse the installed scientific environment before installing another one:

```bash
python -m pip install -e '.[test,campaign,advanced,documents]'
python scripts/run_advanced_campaign.py --stage all
python -m pytest tests/ -v
bash scripts/deploy_phase3.sh
```

The default runner audits retained physical results and executes all 40 solver
fixtures. `--execute-physical --psi4-python /path/to/psi4/python --xtb /path/to/xtb`
explicitly launches sequential native spin, MECP, microsolvation and proton-wire
stages followed by EGNN training. Set `--output /path/to/new-results --scratch
/path/to/new-scratch` for an independent campaign; stage CLI options define budgets
and immutable attempt folders. Failed SCF jobs remain failed records. On Windows,
run `git config core.longpaths true` in this repository before staging the deeply
nested native evidence. This changes only this checkout's Git configuration.
`--kinetics-input complete_model.json` accepts independently certified physical
inputs only after all 20 state/reference and 18 transition-state records pass
the composition, temperature, stationarity, protocol and artifact-hash gates.
Certificates retain an explicit review-attestation boundary: hashes bind bytes,
not chemical truth. Source structures and native calculations remain inspectable.

The spin publication includes inputs, commands, native stdout/Psi4 output,
analytic gradients, state responses, rejected attempts and provenance. Large
numbered Psi4 files `psi.<pid>.97` (DF-SCF B-matrix integrals) and
`psi.<pid>.64` (DIIS storage) stay unchanged on the originating workstation;
they are omitted from Git and the release ZIP, with exact paths, bytes and
SHA256 in `data/phase3/spin/evidence_publication.json`. These are reconstructible
intermediates identified by the [official Psi4 PSIO file map](https://psi4.github.io/psi4docs/master/autodoc_psifiles.html).
Other scientific outputs remain published; the release is not a full scratch backup.

The user-specified baseline is **383.15 K, toluene, 0.05 equivalents tBuOK**;
the base grid is 0.01-0.20 equivalents. Total salt, free tert-butoxide and
neutral-alcohol activity are distinct inputs. The current physical rate network
lacks the required transition-state and off-cycle thermochemistry, so no
computed catalyst activity or lifetime ranking is released.

## Install and verify

Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,campaign]'
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
For the complete campaign tests, set `PINCER_INSTALL_EXTRAS=test,campaign`.
`PINCER_COMMIT_MESSAGE` supplies an explicit campaign message; omitting it keeps
the original core commit convention. Data are included in the staged file set.

See [English report](docs/TECHNICAL_REPORT_EN.md) and
[中文技术报告](docs/TECHNICAL_REPORT_ZH.md) for derivations, citations, algorithms,
evidence limits and the intended connection to pincer catalysis research.

## Actual quantum campaign and evidence

The 24 structures are designed analogues with explicit ancillary ligands and
formal charge/occupation assumptions. Co is modeled as Co(I). The three modeled
activation sites are PNP nitrogen, pyridine-PNN benzylic carbon, and bipyridine
PNN(O) oxygen; a universal N–H activation assignment would be chemically wrong.
GFN2-xTB does not establish their ground spin or experimental accessibility.

The revised conditions are ALPB toluene (`gsolv`), 383.15 K, tBuOK 0.05 equivalent,
with a requested 0.01–0.20 equivalent base grid. Full Hessians provide thirteen
RRHO/qRRHO temperature evaluations. The ALPB electronic potential is held fixed
across these temperatures; solvent temperature derivatives were not calculated.
One ideal 1 atm-to-1 M correction is applied per molecular species at each T.
The contact pair, free anion, and tBuOH are distinct molecular references; their
existence does not determine dissolved salt speciation or free-base activity.

Key source collections:

| Collection | Meaning |
| --- | --- |
| `data/datasets/catalyst_descriptors.csv` | 24 design rows; missing intended structures remain empty |
| `data/structures/` | Portable, exact-byte copies of selected native XYZ geometries and metadata |
| `data/datasets/best_conformers_portable.json` | Electronic-energy selections with original provenance |
| `data/datasets/thermochemistry/` | Fully revalidated certificates, 13T values and balanced reaction free energies |
| `data/campaign/neb/` | Actual paths, force traces, convergence and strict TS acceptance records |
| `data/condensation_search/` | Explicit neutral/anionic cluster attempts and failures |
| `data/ionic_reference_retry/` | Verified K–OtBu contact-pair reference and actual restart audit |
| `data/native_evidence/` | CRC/SHA256-verified ZIP archives of native inputs, outputs, gradients and Hessians |
| `data/datasets/kinetics_readiness.json` | Why the requested real kinetic grid was not executed |
| `data/datasets/surrogate_readiness.json` | Missing accepted barrier-label coverage; no trained model |

The electronic best-found geometry need not be the lowest *certified* free-energy
geometry. The thermal export chooses the lowest certified G at 383.15 K and uses
that same certificate at every temperature; it never attaches its thermochemical
correction to another conformer. Positive curvature alone cannot establish the
intended catalyst identity. In particular, Co bipyridine active rearrangements
are preserved as failures of that intended identity.

Plots in `examples/plots/` have vector SVG/PDF and PNG companions. The Ru/Mn figure
shows certified endpoint thermodynamics only. The requested TOF heatmap filename
contains an explicitly labeled availability plot because no complete kinetic
network is available. Hatched missing cells are not zero activity. There is no
claimed mechanistic rate-controlling step, trained prediction, or metal ranking.

## Replay and new calculations

The actual runtime used native xTB 6.7.1 through an ASE calculator adapter, rather
than an unavailable xtb-python binding. Install a working native xTB separately,
then set `PINCER_XTB` to that executable. Reuse an existing scientific environment
when available. The optional `campaign` dependencies include RDKit, ASE, SciPy,
scikit-learn, and psutil. The `documents` extra supplies the PDF authoring tools.

Revalidate the published evidence and regenerate the figures:

```bash
python scripts/export_thermochemistry_dataset.py
python scripts/export_reference_reaction_thermodynamics.py
python scripts/build_campaign_plots.py
python -m pytest tests/ -v
```

The bilingual PDF renderer uses Windows Times/SimSun/SimHei fonts by default;
provide `--font-directory` with those font files on another host. Rebuild and
render every page for inspection with the `documents` dependencies:

```bash
python scripts/build_monograph_pdfs.py --output-directory work/pdf-rebuild
python scripts/verify_monograph_pdfs.py work/pdf-rebuild work/pdf-qa
```

The second command records both parser page counts, every rendered page hash,
text bounds, and the actual PDF rendering runtime. Inspect the resulting contact
sheets and full pages before claiming visual quality; numeric bounds alone do
not prove legibility or correct layout.

Historical native absolute paths in raw JSON identify the originating workstation.
Derived CSV paths are repository-relative. Certificate exporters resolve the local
same-directory evidence and require its full SHA256. Data have Git text conversion
disabled because changing native line endings would invalidate source hashes.
Each native archive contains work-root-relative member names and an explicit
selection/omission manifest. Binary restart caches, `xtbtopo.mol` topology caches,
`.xtboptok` sentinels and files outside the manifest's explicit selectors are
omitted. Saved XYZ, expected molecular graphs, native WBO, full output, restart
source hashes, commands, final 300 K electronic energies, and force evidence remain.

Start a *new* conformer search under a unique name so published data are preserved:

```bash
python scripts/run_conformer_campaign.py --campaign-name independent_run --hours 3 --workers 2 --threads 2 --solvent toluene --xtb /path/to/xtb
```

The new run writes `data/independent_run/` and separate scratch directories.
`--resume` continues an existing run under its original deadline; it does not
silently allocate a fresh three hours. Select that run's best-index and control
paths. Changing the scientific protocol requires a new campaign name. Select output
paths when launching the Hessian and NEB scripts. Inspect each script's `--help`
for bounded resources and explicit output/scratch paths. Reproduction can change
the best-found conformer with software/CPU differences and does not guarantee TS
convergence. The original failed records remain scientifically meaningful.

The eight-ODE solver and five-fold grouped surrogate are reusable algorithms.
Their analytical regression fixtures are explicitly synthetic software tests,
separate from the real campaign data. Their production interfaces require complete
certified compositions, identities, charges, protocols, free energies, TS modes,
and base-activity assumptions before producing scientific predictions.
