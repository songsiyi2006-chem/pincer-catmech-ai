#!/usr/bin/env bash
# Read-only. Run on the institutional server after legitimate login.
set -u
date -u
uname -a
command -v lscpu >/dev/null && lscpu
command -v free >/dev/null && free -h
df -h .
command -v quota >/dev/null && quota -s
for engine in sbatch squeue sinfo qsub qstat xtb python psi4 pw.x cp2k.psmp vasp_std; do
  command -v "$engine" || true
done
command -v sinfo >/dev/null && sinfo -o '%P %a %l %D %c %m'
# Do not infer job allocation, permitted partition, license or GPU entitlement
# from executable presence. Confirm with the institution before submission.
