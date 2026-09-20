# DFT review resolution / DFT 独立复核闭环

Review date: 2026-09-20. This is a read-only source/record review of the integrated `data/dft_pilot` workflow and `reports/TECHNICAL_REPORT_ZH.md`. Only this resolution document was written. No Psi4 jobs, numerical reruns, or heavy tests were started, and historical review/archived execution snapshots were not edited.

The original findings in [DFT_REVIEW.md](DFT_REVIEW.md) and the wrapper-exit issue identified during this follow-up are **CLOSED** in the final reviewed source. No actionable issue remains open within this review's limited coverage. Closure means the specified code path and recorded validation support the repair; it is not a new independent execution of the test suite or a claim that all possible failures were tested.

## Original findings

| Finding | Final status | Reviewed resolution and evidence |
|---|---|---|
| NaN/nonfinite energies bypass the numeric comparison | CLOSED | The summarizer explicitly rejects nonfinite energy, used HOMO diagnostics and ideal S²; also checks native energy and derived charge difference. `test_nan_energy_cannot_supply_a_label` is recorded as passing. |
| Pair protocol omits Psi4 version and redundant settings/script provenance | CLOSED | Pair validation now checks `psi4_version`; summary checks result settings against `settings.json` and runner config, embedded settings hash, executed child-script hash, and spin reference. Settings/script/version mutation checks are present and recorded as passing. |
| Disappearing child can skip parent kill and lose final record | CLOSED | Child kill exceptions are isolated; parent and `Popen` handle termination are attempted separately. Final waiting/escalation and record writing are in a `finally` path. Mocked disappearing-child and launch-exception tests are recorded as passing. Sampling remains an explicitly non-hard RSS monitor. |
| Unconditional extra preflight event | CLOSED | Extra events are counted only when the external manifest exists, explicitly says not launched, and is not the same preflight already present in runs. The missing-external-file test passes. |
| No actual charge pair can make the sign/unit test vacuously pass | CLOSED | Actual pair arithmetic has its own test and explicitly skips when no pair exists. The archived skip reason states that the physical difference remains unavailable. |
| Duplicate successful method/charge runs silently overwrite the accepted state | CLOSED | The summarizer now raises on an ambiguous duplicate and demands an accepted-run registry. The duplicate-success mutation test passes. |

## New finding

### CLOSED — P2: unsuccessful wrapper runs exited successfully at the shell

The initial follow-up snapshot recorded `timeout`, `memory_limit`, `launch_failed` or `monitor_error` correctly, but reached the final print without returning a failure code. The CLI called `main()` directly, so a handled failure yielded exit status 0. The scientific summary still excluded failed records; this did not change either completed neutral energy, but could mislead a shell-based caller.

Final read-only verification confirms `return 0 if record["status"] == "completed" else 1` and `raise SystemExit(main())`. The existing mocked launch-failure test now asserts a nonzero return code. The refreshed native archive/process-control log records 14 tests, 13 passed and 1 explicitly skipped. Historical executed scripts remain separate from this repaired future-run wrapper. No Psi4 rerun was used for this repair.

## Chinese technical-report cross-check

Within the following inspected records, **no numerical or scientific-boundary discrepancy was found**:

- `data/dft_pilot/summary.json` agrees with two completed neutral runs, two timed-out anion runs, two prelaunch rejections, no complete charge pair, and unavailable method spread. Reported energies are −531.5663974136239 Eh (PBE0) and −532.1735781682701 Eh (B3LYP). Reported run times and RSS values agree to the displayed precision.
- Neutral/anion electron counts and multiplicities are 84/singlet and 85/doublet for the stated C9H8N2O model. The report correctly identifies inherited xTB-generated fixed nuclear coordinates and does not identify them as a confirmed SI substrate, DFT minimum or Cu site.
- `results/neutral_transport_audit.json` agrees with all four grid errors, error ratios, 160-cell normalized error, Radau residual, site/film residuals and forward/reverse integrated-source balances quoted in the report. These remain synthetic neutral software-verification values.
- `results/validation/summary.json` and `tests.txt` report 42 tests, 0 failures, 0 errors and 1 explicit skip: 41 passed plus 1 skipped. The report states this accurately and does not equate the tests with scientific mechanism validation.
- Parsed literature inventory contains 14 source entries, 13 novelty rows and 35 queries. This verifies inventory counts only; the current review did not re-search the web or independently revalidate every literature claim.
- The initial delivery-manifest snapshot listed 58 files and 6,009,327 bytes. After the wrapper/test/log repair, the final reviewed manifest lists 58 files and 6,009,462 bytes; both match the report's approximately 6.01 MB. The scratch manifest lists 8 files and 849,081,395 bytes, matching approximately 849.08 MB. These are decimal MB as written, not MiB. Exact byte totals are snapshot-specific and the manifest identified below is authoritative for this review.
- The report consistently withholds charge-pair energies, reduction potentials, Cu mechanism/ranking conclusions, AI gains, industrial benefits and full TEA/LCA claims. Nonconvergence is not interpreted as a nonexistent anion. Full-cell electrical accounting is separated from FE and the neutral model.

This pass did not repeat the historical xTB audit, rerun environment recovery, independently verify the user-scope history, review every external source, or test every OS/process failure mode. Those subjects are outside the stated limited closure review. No new scientific validation is implied.

## Reviewed source snapshot

| File | SHA256 |
|---|---|
| `data/dft_pilot/summarize_dft.py` | `d00a3d0b92834292fd3a9828dada385b7f334590e1af62732c63386c404763c1` |
| `data/dft_pilot/bounded_runner.py` | `b4891123df75bbe2b212f26cbf5627eab04d8ec7c47ea8a94fc8bf3fff8935c3` |
| `data/dft_pilot/test_dft_archive.py` | `f73049bb164b4bda9341b77139c5c87fb832092b96df1e57511ad554310da8b8` |
| `data/dft_pilot/test_bounded_runner.py` | `904e1d2ad55072229ffddf136cacfcffc668658c75e1b42b3dabeaf3eab08fd1` |
| `data/dft_pilot/delivery_files.json` | `36c297e983a4741ee5bd86c91eb0fcb9ae2e442beef65475e37e26ed884fc70f` |
| `reports/TECHNICAL_REPORT_ZH.md` | `5f29eaf267180898aa52e1da18868655f53681700f61f768cc961a944207f161` |
