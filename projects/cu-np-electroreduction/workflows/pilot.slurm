#!/usr/bin/env bash
#SBATCH --job-name=cu-np-preflight
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=pilot-%j.out
# Template: allocation/partition/module names deliberately not invented.
set -euo pipefail
: "${SLURM_JOB_ID:?Submit through an approved scheduler allocation, not a login-node production run}"
: "${PROJECT_ROOT:?Set absolute path to projects/cu-np-electroreduction}"
: "${PILOT_OUTPUT:?Set a fresh absolute output directory}"
: "${XTB_EXE:?Set verified xtb executable}"
cd "$PROJECT_ROOT"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
python scripts/run_pilot.py --xtb "$XTB_EXE" --threads "$OMP_NUM_THREADS" --output "$PILOT_OUTPUT"
