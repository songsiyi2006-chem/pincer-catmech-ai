# Independent code and scientific-boundary review

Review date: 2026-09-20. Scope: `scripts/science.py`, `run_pilot.py`, `refine_pilot.py`, `summarize_pilot.py`, `process_scenarios.py`, configuration, tests and the 20 archived native calculation records. The reviewer did not rerun quantum calculations. This report does not assess unpublished experimental evidence or unavailable catalyst models.

## Verified current evidence

**EXECUTED RESULT.** The native archive contains 20 unique calculation IDs. Each recorded output checksum matched its retained file; the energies were independently reparsed from the retained stdout/stderr and agreed with the records. All 20 native processes completed; completion does not mean every scientific check passed.

**EXECUTED RESULT.** The original Hessian fails the configured -20 cm^-1 threshold because of its -50.49 cm^-1 mode. That rejected result remains archived. The repaired neutral geometry has 54 vibrational modes with the lowest at 78.64 cm^-1, in addition to six projected zero modes. Its nuclear composition is C9H8N2O (20 atoms); neutral and anion single points use identical coordinates and correct 84/85-electron parity.

**EXECUTED RESULT.** Independent recomputation of paired charging energies reproduces the summary. The GFN1/GFN2 differences are about 0.782 eV in gas and 0.832 eV in ALPB(DMF), both above the declared 0.2 eV interpretation stop threshold.

**MODEL INFERENCE.** The interpretation is appropriately restricted: equilibrium ALPB fixed-nuclear-geometry charging energies are not nonequilibrium electron affinities, calibrated electrode potentials, transfer barriers or evidence for a Cu/P catalytic effect. Neither a catalyst atom nor an explicit electrochemical interface appears in this pilot.

## Required regeneration repairs identified

1. **Failure-aware summary generation.** At the reviewed revision, `summarize_pilot.py` reparses every record without branching on recorded status or catching parse errors. A retained timeout/nonconverged calculation therefore stops the entire summary before its failure can be indexed. Preserve a row and raw-file index for every record; set unavailable energies to null, retain failure reasons, and restrict pair calculations to completed, successfully parsed jobs. Hash mismatch should remain a hard error. This defect does not alter the present energies because all 20 native calls completed.

2. **Derive gate state from evidence and configuration.** The reviewed summary hardcodes `initial_minimum_rejected=True` and 0.2 eV instead of deriving them from the actual Hessian and `config/pilot.json`. It also does not explicitly require the repaired Hessian to pass before producing derived charging comparisons. Read the threshold from configuration and make the revised-minimum check a derivation gate. This protects future regeneration; the present raw Hessians and configured 0.2 eV do support the current flags.

These recommendations were sent to the lead for implementation. See the review addendum below for the final status after changes.

## Other bounded limitations

- **MODEL INFERENCE:** Runner and summarizer use different Hessian representations: 54 symmetry-labelled vibrations versus all 60 modes. Both correctly classify the present structures; the revised summarizer now reports both counts and the 78.64 cm^-1 vibrational minimum separately.
- **MODEL INFERENCE:** The execution wrapper preserves failed/running directories and will not overwrite them; it is conservative and idempotent for completed archived jobs. It does not automatically retry a failed ID. Recovery requires a documented new batch/ID. A directory created before its first `record.json` write would need explicit quarantine/recovery.
- **MODEL INFERENCE:** Cached generated geometry is checked by checksum, but changing the configured scaffold to another isomer with the same formula could reuse the old cached structure unless its source SMILES/seed/provenance are also checked. This does not affect the fixed 1a pilot. Future configurable use should reject such a provenance mismatch.
- **MODEL INFERENCE:** Methyl-rotor repair uses fixed atom indices and assumes this exact 20-atom ordering. It should remain a named one-off reconstruction script, not a general conformer-search workflow. A mode eigenvector analysis was not independently demonstrated here; the successful new Hessian validates the resulting local minimum, not a mechanistic assignment to the original mode.
- **MODEL INFERENCE:** Several numerical utilities do not reject all NaN/infinity inputs, and the ideal-site network uses directly scaled tree weights. These are not used with pathological values in current results. Production fitting would need finite-input checks and numerical scaling over its admitted rate range.
- **MODEL INFERENCE:** The grouped-split helper prevents one identity being assigned to inconsistent groups. It cannot itself prove that a user-defined identity corresponds to an independent chemical scaffold. The current prospective file reserves three condition groups; there is no fitted ML model or prediction-validation claim.

## Formula and conservation checks

**EXECUTED RESULT.** Additional tests use independent scientific relationships rather than matching output strings:

- Hartree-to-eV-to-kJ/mol conversion and paired energy subtraction.
- Explicit element/charge balancing with electrons and protons.
- Potential-reference round trip and a shared-reference equality; no uncalibrated DMF-to-aqueous conversion was performed.
- Shared-transition-state rate pairs reproduce detailed balance.
- A nondegenerate three-state network reproduces Boltzmann equilibrium, and a driven network gives identical flux on its three edges.
- Full-cell electrical energy computed from an explicit one-kilogram product/electron balance agrees with `U*z*F/(3600*M*FE*recovery)` for M in g/mol.
- The synthetic energy scenario keeps z and product molecular weight unverified; it does not fabricate a measured kWh/kg value.

**EXECUTED RESULT.** Test command: `python -m unittest discover -s tests -v`. At the first independent review, 18 tests passed. After failure-injection regressions, lead repairs and the separately authored kinetic-diagnostics tests were added, the full suite finished with **31 tests passed** in the independent rerun. Tests of mathematical relations are labelled synthetic software tests; archival tests read actual retained results. Passing these tests does not establish electrochemical accuracy, reaction mechanism, selectivity or experimental validation.

## 中文结论

**EXECUTED RESULT：** 实际20次量化调用的原始文件、哈希、能量和电子数通过独立审查；初始负频失败得到保留，修复后最低振动频率78.64 cm^-1。两种半经验方法对充电能的差异约0.8 eV，触发预设解释停止条件。

**MODEL INFERENCE：** 当前结果只能支持“分子预检存在明显方法敏感性”，不能支持“P提升Cu催化选择性”或任何定量工业结论。汇总器的失败任务保留、Hessian成功门槛与配置阈值整改已复核，最终31项测试通过。测试通过与科学机理得到验证是不同层次的证据。

## Review addendum

**EXECUTED RESULT.** Both required findings were repaired and independently reinspected. Failed records remain indexed with null energies; failed startup with no logs remains visible; missing logs on a purportedly completed calculation reject its parsing. Missing files that were already hashed in an archive remain an integrity error. Hessian-derived comparisons now require completed native status, a parsed energy and a passing spectrum. The Hessian and method-spread thresholds are read from configuration, and the initial rejection is derived from the actual archive.

**EXECUTED RESULT.** Four explicitly synthetic failure/configuration tests run on temporary copies: native failure with a misleading energy-looking line; startup failure before log creation; a failed Hessian that nevertheless left a positive spectrum; and a changed method-spread threshold. They all pass after repair. The third test initially failed, exposing a second-order gate bug that the lead corrected. Production raw outputs were not modified by these tests, and no quantum jobs were rerun.

**EXECUTED RESULT.** A scientific-boundary scan of the English technical report confirmed that catalytic predictions remain unclaimed. A scenario-count error was corrected from 36 to 18 (3 voltage values x 3 FE values x 2 recovery values). Remaining limitations above concern future generalization and do not invalidate the stated 20-call molecular preflight.
