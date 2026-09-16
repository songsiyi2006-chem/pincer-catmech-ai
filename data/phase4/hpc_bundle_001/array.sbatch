#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192M
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
set -euo pipefail
: "${HPC_PYTHON:?set absolute site Python 3.10+ path}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec "$HPC_PYTHON" bridge.py run-task --bundle . --index "${SLURM_ARRAY_TASK_ID:?}"
