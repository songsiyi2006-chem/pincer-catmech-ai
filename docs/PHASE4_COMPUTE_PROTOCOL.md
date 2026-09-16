# Phase 4 local computation: diagnosis before expansion

This protocol is a bounded feasibility and numerical-recovery pilot. It is not a production DFT study, a discovery of spin crossover, or a statement that a Nature paper is ready. The relevant historical evidence stays unchanged under `data/phase3/`.

## What failed, and what that does and does not establish

The completed Phase 3 vertical screen scheduled 54 states across 18 Fe/Co/Mn designs. Six source structures were unavailable at the state-slot level; 48 state calculations were attempted and eight reached the stated SCF/gradient finish condition. A timeout is a resource outcome, not proof that an electronic state does not exist.

The Fe bipyridine PNN-OH iPr active geometry is particularly useful for diagnosing the local bottleneck. Its original PBE/STO-3G singlet completed. The triplet completed too, but its reported expectation value was <S²> = 2.6614457635, versus the nominal triplet value 2. Therefore it failed the unchanged absolute spin-contamination screen of 0.1. Allowing more wall time alone does not resolve a converged contaminated solution.

The completed bounded quintet attempt reached SCF iteration 26 before its 120 s timeout. Its last printed energy change was −3.23503 × 10⁻⁸ Eh and commutator RMS residual 5.19880 × 10⁻⁷. The required energy and density/residual tolerances remain 10⁻⁸ and 10⁻⁶, respectively. This was close enough to motivate one longer run, while leaving open whether its final gradient and spin quality will pass. An earlier interrupted attempt of this state stopped near the SAD initialization; it must not be confused with the later 26-iteration record.

The original def2-SVP singlet resource pilot was stopped after 90 s and only five printed SCF iterations. Its energies were still fluctuating. The present evidence supports “insufficient bounded run”, rather than “def2-SVP cannot converge”. A larger basis, by itself, is not proof of method reliability.

The neutral 43-atom hemiaminal + one tBuOH path had two initial 100-step CI-NEB attempts and one 200-step continuation. None met the 0.07 eV Å⁻¹ band-force requirement. The continuation ended at 1.206630 eV Å⁻¹ and showed appreciable oscillation over its last steps. Its sampled maximum electronic energy was not a stationary transition-state energy. Adding the same number of steps again is not the highest-information local experiment. The next path study should first inspect image spacing, proton placement, the reaction coordinate and competing stepwise versus concerted pathways. Any dimer search launched from an unconverged band would be a new independent saddle search; it cannot inherit a “converged NEB” label.

## Fixed pilot, selected before the new results

The ordered structure is `data/structures/Fe_bipyridine_pnnoh_iPr/active/best_found.xyz`; metadata and source bytes are captured before launch. This is a generated and xTB-relaxed model, not an experimentally determined complex. Its metadata identifies the responsive site as hydroxypyridine O-H/O⁻. The earlier generic NH-deprotonation hypothesis must not be silently applied to this scaffold. Chemical identity, proton location, charge balance and the intended experimental structure must be reviewed separately for each ligand family.

| Job | Geometry | Electronic method | Multiplicity | Maximum worker time |
| --- | --- | --- | --- | --- |
| `sto3g_quintet_recovery` | Unchanged mapped active geometry | PBE/STO-3G | 5 | 240 s |
| `def2svp_singlet` | Same geometry | PBE/def2-SVP | 1 | 450 s |
| `def2svp_quintet` | Same geometry | PBE/def2-SVP | 5 | 450 s |

All three use gas phase, density fitting, SAD initialization, 35 radial × 110 spherical DFT points, late SOSCF switching at 10⁻⁴, at most 120 SCF iterations, E convergence 10⁻⁸ and D convergence 10⁻⁶. SOSCF inherits five maximum microiterations from the existing provider. Singlet uses RKS and quintet uses UKS. The derivative is analytic. The very coarse integration grid and the basis choices are deliberate feasibility settings; neither is a production recommendation.

One worker runs at a time. Each worker uses two threads and 500 MiB configured Psi4 memory. The external process monitor enforces 1536 MiB resident memory, including child processes if present, and a host free-memory reserve. The quantum worker deadline is at most 1200 s for the entire run, with the per-job ceilings above. Evidence finalization and hashing may require a small amount of additional time; the receipt records actual elapsed time. Jobs are not submitted if insufficient time or memory remains. A killed worker's partial log remains evidence of failure, not a source of an accepted energy.

The reaction baseline of 383.15 K remains project context only. These are vertical electronic energies, with no thermal corrections, standard-state transformation, solution dielectric or explicit tBuOK. Calling them “383 K solution free energies” would be incorrect.

## Execution and evidence

Use the already installed runtime, from the repository root:

```powershell
& '..\run_campaign.ps1' -PythonArguments @(
  'scripts/run_phase4_local_pilot.py',
  '--python','C:/Users/HUIWEI/miniconda3/envs/phase7/python.exe',
  '--wall-seconds','1200'
)
```

The default output is `data/phase4/local_pilot_001`; scratch is the sibling workspace `../phase4/local_pilot_001`. Both must be new and non-nested. Existing results are never overwritten, and Phase 3 is explicitly excluded as a destination. A repeat requires two different fresh paths:

```powershell
& '..\run_campaign.ps1' -PythonArguments @(
  'scripts/run_phase4_local_pilot.py',
  '--python','C:/Users/HUIWEI/miniconda3/envs/phase7/python.exe',
  '--output','data/phase4/local_pilot_002',
  '--scratch','../phase4/local_pilot_002'
)
```

The run snapshots its driver, the unchanged Phase 3 provider, input XYZ and metadata. It hashes the interpreter, Psi4 package manifest and core binary when present. Every job saves exact JSON input, executable argument vector, combined worker stdout/stderr, native Psi4 output, response and acceptance diagnostics. Native scientific text and arrays are copied into the publication output with SHA-256 verification. PSIO binary intermediates remain at their local scratch paths, indexed with their own byte counts and hashes; they are not silently deleted or represented as included publication files. Output and scratch paths in native records describe the actual historical run and may require relocation after downloading a release.

The prospective unit tests exercise immutable path handling, geometry and protocol mismatches, nonfinite results, spin-contamination rejection and evidence-copy integrity. Their numbers are software fixtures only. Passing them does not validate any electronic-structure approximation.

## Decision after the pilot

1. A completed quintet with acceptable <S²> provides a recovered vertical diagnostic at the stated coarse protocol. It does not establish the ground state, orbital stability, multireference character, a stationary catalyst minimum, or an MECP.
2. A matched pair of converged def2-SVP calculations allows a same-geometry electronic difference to be inspected. It still is not a validated production spin gap. Comparisons must retain the same method/basis/grid/geometry and state identity; no mixed-level subtraction is permitted.
3. A large or sign-changing difference upon protocol changes is an uncertainty warning to investigate, not a reason to select whichever protocol gives the preferred mechanism.
4. Timeouts, memory stops and spin-contaminated states remain failed or ineligible. No threshold is relaxed to increase the success count. The local driver sets `accepted_chemical_label=false`, `production_spin_gap_validated=false`, and `accepted_MECP_count=0` regardless of numerical success.

## What must move to professional electronic-structure/HPC work

First review the exact proposed catalyst, coordination sphere, atom mapping, charge, electron count and proton-responsive site with the experimental chemist. Then select a limited set of chemically motivated competing electronic solutions, including broken-symmetry alternatives where appropriate, and check their orbital stability and identity. Geometry optimization and vibrational characterization are needed before interpreting state minima. A systematic basis/grid/functional sensitivity assessment should be chosen for this transition-metal chemistry; model agreement is stronger when independently justified, rather than obtained by trying many methods until one supports a preferred answer.

For the physical solution mechanism, retain the user conditions: toluene, 383.15 K and total tBuOK loading 0.05 equivalents as a baseline. Determine whether potassium coordination, ion pairing, aggregation, residual water and tBuOH activity affect the key equilibria. Total salt loading is not automatically the activity of free tBuO⁻. Standard-state and thermal corrections must use one documented convention.

For any proposed transition state, require the appropriate stationary force tolerance, a consistent Hessian and one reaction-relevant imaginary mode; confirm connections using IRC or a justified alternative and verify the resulting basins. A spin-crossing hypothesis additionally needs same-Hamiltonian, same-geometry state energies/gradients, an actual constrained crossing search and investigation of whether intersystem crossing can compete kinetically. A crossing by itself does not determine the transition probability or catalytic rate.

Only after a credible, mass-balanced free-energy network exists should physical microkinetics and uncertainty propagation be presented as chemical predictions. Preserve separation between that future network and the current software-fixture kinetic tests. Publish the relevant negative results and limitations rather than assigning a journal level to calculation counts.
