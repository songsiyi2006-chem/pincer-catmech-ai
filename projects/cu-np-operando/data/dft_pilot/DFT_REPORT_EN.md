# Local DFT diagnostic report

**Evidence class: CALCULATION.** 2 completed; 2 failed/terminated; 2 prelaunch resource rejections.

This is a bounded, real calculation on independently generated C9H8N2O. It is not a completed CuN4/CuN3P1 catalytic study. The accepted prior GFN2-xTB neutral geometry was used as fixed input; no original SI coordinates or Cu active sites are implied.

| Run | Method/charge | Status | Energy / Eh | Wall / s | Peak RSS / MiB |
|---|---|---|---:|---:|---:|
| DFT001-pbe0-neutral | pbe0, q=0 | completed | -531.566397413624 | 108.66 | 612.52 |
| DFT002-pbe0-anion | pbe0, q=-1 | timeout | unavailable | 300.06 | 622.55 |
| DFT003-b3lyp-neutral | b3lyp, q=0 | resource_preflight_not_launched | unavailable | not launched | 0.00 |
| DFT003b-b3lyp-neutral | b3lyp, q=0 | completed | -532.173578168270 | 117.23 | 612.55 |
| DFT004-b3lyp-anion | b3lyp, q=-1 | timeout | unavailable | 300.13 | 622.18 |

## Charge-pair result

No complete charge pair: ΔE, attachment energy and inter-functional spread are **unavailable**.

Only a completed neutral/anion pair with identical nuclear hashes, basis content, Hamiltonian, DF settings, grid and convergence tolerances supplies a difference. Restricted/unrestricted spin treatments differ as required by charge state. A documented SOSCF rescue, if present, changes the iterative solver only. No last SCF iteration from a timeout is accepted as an energy label. No xTB total energy appears in the DFT differences.

## Interpretation and limits

At fixed nuclei R, ΔE = E(anion; R) − E(neutral; R), and the electronic attachment energy is A = −ΔE. Positive A means a lower anion energy within the finite-basis model; negative A means the charged state is higher. Neither establishes a calibrated physical electron affinity or an electrode potential.

SCF convergence is not SCF stability. No orbital-Hessian stability analysis was performed, and a single SAD initial guess was used for the standard runs. A second-order SCF solver does not substitute for a stability analysis. Spin expectations, if a run completes, describe the corresponding Kohn–Sham determinant. Positive anion occupied orbital energies are warnings; finite diffuse basis functions can confine a state that is not physically bound. A negative orbital energy would likewise not prove binding.

The single def2-SVPD basis and 75×302 integration grid were not convergence-tested. No DFT optimization, vibrational Hessian, non-equilibrium solvent polarization, electrode electron reservoir, potential-reference calibration, Cu coordination site, kinetic network or experiment was included. Agreement between two approximate functionals would not eliminate their shared systematic error.

## Resource and next-decision evidence

The native inputs, output, software version, spin populations, final energies, orbital coefficients (where completed), wall time, sampled memory history and hashes are preserved. Each quantum launch used two threads, a 512 MiB Psi4 memory setting, and a 300 s timeout. The process-tree RSS cap was tightened as host memory changed; no quantum jobs ran in parallel. Prelaunch RAM failures are recorded separately from native runs. A 0.1 s RSS monitor is a sampled guard, not an OS-enforced absolute memory quota.

The next defensible work is a separately budgeted SCF/basis/anion-binding diagnostic or a sourced electrode model after the missing SI/HPC requirements are resolved. These results do not support catalyst ranking, reduction potentials, selectivity or industrial productivity. See [Psi4 SCF/stability documentation](https://psicode.org/psi4manual/master/scf.html#stability-analysis) for the distinction between a converged SCF solution and a stable one.

## Reproduce

Run `python summarize_dft.py`, followed by `python -m unittest discover -s . -p "test_dft_archive.py" -v`. Rebuild the bilingual notes with `python write_summary_reports.py`. These commands require no Psi4 calculation. The computational rerun command and exclusion policy are in README.md. Rebuildable scratch is retained locally and listed by hash, but excluded from Git.
