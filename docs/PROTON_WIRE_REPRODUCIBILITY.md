# Proton-wire execution and evidence finalization

The configurable driver and published finalizer are a prospective improvement made after the archived physical campaign. No quantum calculation was repeated while implementing these CLI and packaging changes. The existing campaign's authoritative execution snapshots are:

| Executed stage | Archived source | SHA256 |
| --- | --- | --- |
| Two original seed runs | `data/phase3/proton_wire/driver_initial_run.py` | `67afe36bfd92909da5d66b9336f07ae977c00956eb0000acbff10b2676874bcb` |
| Immutable continuation | `data/phase3/proton_wire/driver_continuation_run.py` | `f37a25cc24239a790274d9290be11da54e56902b8e146e7f4952b3f183fe6b0a` |

The current source in `scripts/run_phase3_proton_wire.py` must not be described as the source that produced those archived quantum results. Their ZIPs, summaries and verification receipt remain unchanged. The existing receipt SHA256 is `2195539f145aa84e54531599857be4d86e24d2eda3087b510ec59dc78b9c76d6`.

The following commands are a recipe for a future explicitly requested physical run, not commands executed during this reproducibility update. Run them from the repository root after activating the existing campaign Python environment with this package installed. Choose fresh, separate, non-nested output and scratch directories. The cluster source may be the existing mapped Phase3 clusters or a fresh master-campaign solvation output. The xTB path below is the originating workstation's installed executable; substitute your existing executable on another computer. These commands use the selected Python directly and do not depend on an unpublished wrapper.

```powershell
$wireOutput = '..\phase3_rerun\proton_wire'
$wireScratch = '..\phase3_rerun_scratch\proton_wire'
$wireClusters = 'data\phase3\solvation\clusters'
$wireXtb = 'C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin\xtb.exe'
$wirePython = (Get-Command python -ErrorAction Stop).Source
$wireLibrary = Split-Path -Parent $wireXtb
$env:PATH = "$wireLibrary;$env:PATH"
$env:OMP_NUM_THREADS = '2'
$env:MKL_NUM_THREADS = '2'
$env:OPENBLAS_NUM_THREADS = '2'
$env:PYTHONUTF8 = '1'

# Two neutral 43-atom seeds, 100 FIRE/NEB steps each, one shared wall deadline.
& $wirePython scripts/run_phase3_proton_wire.py --xtb $wireXtb --output $wireOutput --scratch $wireScratch --cluster-source $wireClusters --budget-seconds 1200 --neb-steps 100

# Optional: one immutable continuation, before finalization and while time remains.
# This reads and preserves the first run's original deadline; it does not reset it.
& $wirePython scripts/run_phase3_proton_wire.py --xtb $wireXtb --output $wireOutput --scratch $wireScratch --continue-first --neb-steps 200

# No quantum jobs: archive every raw file, verify hashes, and produce CSVs/receipt.
& $wirePython scripts/run_phase3_proton_wire.py --output $wireOutput --scratch $wireScratch --finalize-evidence

# Read-only archive audit; an existing verification receipt is not rewritten.
& $wirePython scripts/run_phase3_proton_wire.py --output $wireOutput --scratch $wireScratch --verify-evidence
```

Check each physical command's summary before choosing the optional continuation. A failed search is preserved as a scientific negative result; it is not converted into a software exception or a successful TS. A continuation is refused after finalization, after a previous continuation, after an accepted first saddle, or when less than 60 seconds remain in the original allowance. Finalization remains available when the quantum deadline has elapsed, because it performs only local file inspection and compression. `--continue-first`, `--finalize-evidence` and `--verify-evidence` are separate mutually exclusive modes.

The standalone published helper provides the same finalization operation:

```powershell
& $wirePython scripts/finalize_phase3_proton_wire.py --output $wireOutput --scratch $wireScratch
& $wirePython scripts/finalize_phase3_proton_wire.py --output $wireOutput --verify-only
```

Prospective physical runs snapshot their exact initial/continuation driver source before computation in `run_sources/`, with hashes in `summary.json`. Finalization verifies those snapshots and records its own source separately. It never assumes that the currently installed driver produced an older calculation.

The finalizer archives every scratch file, including electronic restart/cache files, into immutable numbered ZIPs grouped by path run and split below 95 MiB. It records each archive name, member path, byte count and SHA256 in one combined manifest. It verifies archive membership, checksums, endpoint graph/Hessian evidence and exact restart coordinates, then exports the sampled-band and endpoint-thermochemistry tables. Incomplete or failed paths remain represented. A passed finalization receipt means the evidence checks passed; it does not mean a TS was accepted.

Reinvoking finalization on an already passed receipt performs a read-only archive audit and leaves the receipt and summaries unchanged. Archive corruption is an error. No raw file is deleted. A failed packaging attempt retains its staging directory for inspection; an overlarge individual compressed file produces an explicit error rather than dropping data.

Validation of this prospective tooling uses synthetic byte fixtures only for archive splitting, corruption detection, source-provenance boundaries, read-only repeat finalization, output/scratch separation and custom-directory CLI execution. The four existing atom-mapping/model-rejection tests remain separate. These tests do not fabricate electronic energies, gradients, Hessians or saddle evidence.
