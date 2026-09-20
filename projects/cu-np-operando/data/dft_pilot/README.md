# Bounded molecular DFT diagnostic

This directory is a **real local calculation archive**, not a Cu-catalyst or electrochemical prediction. The molecule is independently generated 3-methylquinazolin-4(3H)-one (C9H8N2O), inherited from the prior molecular xTB preflight. Its nuclear coordinates are fixed throughout all new calculations.

The intended matrix is PBE0 and B3LYP, each with neutral singlet and anion doublet, at def2-SVPD. A job counts as completed only when the native SCF log, structured result and resource-run record agree. Incomplete output does not supply an energy label.

## Reproduce without running quantum chemistry

Use ordinary Python with its standard library for the archive checks; the separate runner-control tests also use the existing `psutil` dependency:

```bash
python summarize_dft.py
python -m unittest discover -s . -p "test_dft_archive.py" -v
python -m unittest discover -s . -p "test_bounded_runner.py" -v
```

These commands reparse the archive, verify SHA256 hashes, confirm native convergence/final energy, validate electron/spin counts, check identical geometry/basis/SCF protocol for each neutral/anion pair, and regenerate `summary.json`. The tests mutate temporary copies to verify fail-closed behavior. They do not modify native calculations.

The final run status is **two completed neutral calculations and two timed-out anion calculations**, plus two prelaunch memory rejections. Consequently no complete charge pair or attachment energy is available. The charge-pair integration test is explicitly skipped for this reason. No SOSCF rescue was run; its optional runner flag is an unexecuted future capability.

The final workflow code includes repairs from independent review: finite-value/protocol/hash guards, rejection of ambiguous duplicate successful states, event-derived resource counts, and reliable per-process termination/final records. Historical `executed_*.py` snapshots remain unchanged and identify the actual code used for each calculation.

## Run a new bounded quantum job

The runner interpreter requires `psutil`; the child interpreter requires Psi4. Reuse an existing environment. The new output directory must not exist:

```bash
python bounded_runner.py \
  --psi4-python /path/to/psi4/python \
  --geometry runs/DFT001-pbe0-neutral/input.xyz \
  --method pbe0 --charge 0 --rss-limit-mib 896 \
  --output runs/NEW-pbe0-neutral
```

For the charged state, use `--charge -1`; for the second functional, use `--method b3lyp`. The job uses two threads, a 512 MiB Psi4 internal allocation request and a 300 s wall-time limit. The sampled process-tree RSS cap is configurable only to 768/896/1024 MiB; the runner requires twice that cap in available RAM before launch. The cap is polled every 0.1 s, so it is **not** a kernel-enforced allocation guarantee and transient spikes could be missed. Jobs were run sequentially.

`memory_preflight_rejection.json` preserves a prelaunch refusal when available memory fell below twice a 1 GiB cap. That event is not a quantum calculation. Subsequent jobs can use a tighter cap after the measured neutral benchmark demonstrated sufficient room.

## Scientific definition and limits

At the identical fixed nuclear geometry R, for one method and basis:

`delta_E = E(anion; R) - E(neutral; R)`

`electronic_attachment_energy = E(neutral; R) - E(anion; R) = -delta_E`

Positive attachment energy means the anion is lower **within that finite-basis electronic model**. Negative attachment energy means the computed charged state is higher than the neutral state at the same nuclei. Neither quantity is an electrode reduction potential, a reaction barrier or a free energy. No xTB energy is combined with a DFT energy in these differences.

def2-SVPD includes diffuse augmentation, but one diffuse basis does not establish anion binding. Positive occupied Kohn–Sham orbital energy is an additional warning, and even a negative value would not by itself prove a stable physical gas-phase anion. Electron leakage, basis confinement, self-interaction/delocalization error and SCF solution dependence need further examination before a quantitative electron-affinity claim. The archived orbital energies and Molden file enable that follow-up; no calibrated bound-state classification is supplied here.

SCF convergence means a self-consistent solution was reached, not that the electronic solution is stable. This pilot did **not** perform an orbital-Hessian stability analysis or compare multiple SCF initial guesses. No DFT geometry optimization, vibrational Hessian, basis/grid convergence, solvent calculation, electrode/reference calibration, Cu active-site calculation or experiment was performed. Two approximate functionals may share systematic errors; agreement is not calibration.

Psi4's [SCF documentation](https://psicode.org/psi4manual/master/scf.html#stability-analysis) distinguishes convergence from stability, and its [basis tables](https://psicode.org/psi4manual/master/basissets_tables.html) and [functional table](https://psicode.org/psi4manual/master/dft_byfunctional.html) identify the implemented method families. Reproduction should preserve the archived Psi4 version and native functional definition rather than assume all programs use identical B3LYP conventions.

## Archive policy

Every run retains `input.xyz`, settings, the executed child/runner source snapshots, native `psi4.out`, stdout/stderr, resource history and hashes. Successful runs additionally retain machine-readable energy/spin/orbital results and a Molden file. Failure/timeout records are retained and are never used as training labels.

`excluded_scratch_manifest.json` records relative paths, sizes and hashes of any remaining `scratch/` or root-level `psi.*` rebuildable intermediates. Keep them in the local work area; do not add them to Git. Their absence in a redistributed archive does not invalidate the preserved input/native-output evidence. No scientific result depends on those intermediates.
