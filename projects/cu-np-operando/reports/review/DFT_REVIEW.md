# Independent read-only DFT workflow review

Reviewed on 2026-09-20. Scope: `work/compute_audit/dft_pilot/psi4_single_point.py`, `bounded_runner.py`, `summarize_dft.py`, `test_dft_archive.py`, and current provenance/summary records. No Psi4 jobs or heavyweight tests were run. Findings below refer to the source snapshot read in this review, and may be resolved by the owning agent subsequently. Historical `executed_*.py` files and their hashes must remain unchanged when current workflow source is repaired.

## Assessment

The existing path correctly separates resource rejection, failed/timed-out calculation, and completed calculation. `fail_on_maxiter=True` prevents a routine SCF iteration-limit failure from intentionally returning a usable result. The summary requires completed runner status, completed child status, a native SCF convergence message, and a matching native final energy. A timeout archive is excluded even if partial SCF output exists. The molecular electron count (C9H8N2O: 84 neutral / 85 anion), charge/multiplicity convention (0/1 and -1/2), and Hartree-to-eV/sign conventions are internally consistent. The scope correctly avoids reduction-potential, stationary-point, physically bound-anion, or Cu-interface claims.

The following fixes are recommended before declaring the validator robust. These are implementation findings, not claims that the recorded finite completed energy is wrong.

## Actionable findings

### P2 — Nonfinite energies can pass validation

Location: `summarize_dft.py:79`, energy comparison and later pair subtraction.

`abs(native_energy - result_energy) > tolerance` is false when `result_energy` is NaN. Python's JSON parser accepts NaN by default. A corrupted/nonfinite result can therefore pass the comparison and propagate a NaN label or method spread. `s2_orbital_overlap` is checked for finiteness, but energy and reported orbital diagnostics are not.

Fix: require `math.isfinite` for the parsed native energy, result energy, used orbital diagnostics and final derived differences; reject nonfinite JSON values. Add an archive-mutation test that changes energy to NaN, refreshes only the enclosing result-file digest, and confirms rejection before any label is emitted. No electronic-structure rerun is needed.

### P2 — Pair provenance does not enforce the full recorded protocol

Locations: `summarize_dft.py:19–35`, `summarize_dft.py:81`.

The pair guard checks method, basis, solvent, geometry flags, identical coordinate hash, most SCF options and orbital-basis file hash. This is useful. However, it omits `psi4_version`. Different Psi4 versions can be treated as the same protocol even if implementation/default details differ. The summary also uses `result["settings"]` without checking it equals the separately archived `settings.json` and runner `config`, and it does not verify the embedded `settings_sha256` or `child_script_sha256`. File-level hashes alone do not prove those mutually redundant fields agree.

Fix: verify embedded settings/script hashes against the immutable archived copies; compare parsed settings across result, settings file and runner record; compare Psi4 version across paired calculations. Explicitly verify the expected RKS/UKS reference for the stated multiplicity while allowing the documented SOSCF iterative rescue. Add version/settings mismatch tests. Preserve existing executed source copies rather than editing them retroactively.

### P2 — Termination can skip the parent kill and lose the final record

Location: `bounded_runner.py:89–99`.

One `try` block encloses all child kills and the parent kill. If a child exits between enumeration and `child.kill()`, `NoSuchProcess` jumps past `parent.kill()`. The subsequent `proc.wait(timeout=10)` may then raise `TimeoutExpired`, so the requested timeout/memory termination is not reliably enforced and the final `run_record.json` may never be written. Launch/monitor errors also do not have a final-record `finally` path.

Fix: catch disappearance per child; always attempt the parent termination separately, retry/escalate after bounded waiting if appropriate, and write failure/termination metadata in a `finally` path. Record monitor/launch exceptions instead of turning them into success. A mocked disappearing-child test can verify this without executing Psi4. The 0.1 s RSS sampling interval remains a sampled limit, correctly not a hard allocation quota.

### P3 — Preflight-rejection count includes an unconditional extra event

Location: `summarize_dft.py:129`.

`1 + sum(...)` assumes an external preliminary rejection always exists. The test archive copies only `runs/`; the extra file is absent there, yet the function still invents one preflight event. It also makes the summary nonportable if an archive is delivered without that file.

Fix: load `memory_preflight_rejection.json` only when it exists, confirm its explicit unlaunched status, and count it once. Avoid double counting an event already represented in `runs/`. A missing-file test should produce zero extra rejections.

## Test-coverage caveat

`test_native_hash_energy_electron_and_protocol_validation` iterates over `energy_pairs` without requiring that a pair exists. In the summary snapshot initially read here, there were no completed pairs, so its sign/unit portion was vacuous. This is appropriate for an incomplete live archive only if reported explicitly. Once the intended charge pair has finished, test its actual pair arithmetic. Otherwise mark the integration test skipped with a clear reason and separately test the unit-conversion arithmetic using a clearly labelled numerical validation fixture. Do not create a fake physical pair to make the test pass.

The current dictionary `validated[(method, charge)]` silently keeps the lexicographically last successful run if duplicates arise. Current scope may have only one accepted run per key. Before future expansion, define an explicit accepted-run registry or reject ambiguous duplicate successes; iteration order should not become the scientific selection policy.

## Scientific limits retained

- Same-nuclei gas-phase energies at an inherited xTB-generated neutral geometry are electronic attachment diagnostics, not relaxed electron affinities or reduction potentials.
- A converged SCF and reasonable S² do not establish wavefunction stability, a physical bound anion, or basis/grid convergence.
- SOSCF can be an iterative rescue at unchanged Hamiltonian/tolerances; it must remain visible in the record and does not substitute for checking distinct electronic solutions.
- Current methods do not certify Cu coordination, electrode potential, solvent/ion effects, transition states, reaction pathways, kinetics, rankings or industrial performance.

No further chemistry calculation is requested by this review. The identified repairs and tests are archival validation and process-control work.
