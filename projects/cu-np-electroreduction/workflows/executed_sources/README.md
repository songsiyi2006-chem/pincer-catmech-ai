# Exact executed source snapshots

These byte snapshots match `runner_sha256` and `parser_sha256` in the two
native batches. They preserve the original implementation, including the
initial Hessian count check and the then-unused NumPy three-state solver.
They were recovered from the recorded edits and accepted only after exact
SHA256 equality to the per-calculation records. Current corrected code is
under `scripts/`. These copies are provenance, not recommended executables;
their relative-root logic assumes placement as the original scripts.

The refinement driver additionally records its seed transformation in every
R000 input provenance (methyl hydrogens rotated +60 degrees). Raw coordinates
and settings, which fully define every native calculation, are retained.
