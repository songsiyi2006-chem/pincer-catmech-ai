#!/usr/bin/env bash
# Rebenchmark the entire numerical grid and validate native receipts before
# the existing clean-main, full-pytest, conventional-commit and no-force push.
set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd -- "$root"
# Perform the upstream/branch preflight before regeneration dirties the tree.
command -v git >/dev/null || { printf 'ERROR: git is required.\n' >&2; exit 1; }
[[ "$(git symbolic-ref --quiet --short HEAD)" == main ]] || { printf 'ERROR: Checkout main first.\n' >&2; exit 1; }
[[ ! -e "$(git rev-parse --git-path MERGE_HEAD)" && ! -d "$(git rev-parse --git-path rebase-merge)" && ! -d "$(git rev-parse --git-path rebase-apply)" ]] || { printf 'ERROR: Finish the current Git operation first.\n' >&2; exit 1; }
git fetch --prune origin
if git show-ref --verify --quiet refs/remotes/origin/main && ! git merge-base --is-ancestor origin/main HEAD; then
    [[ -z "$(git status --porcelain)" ]] || { printf 'ERROR: Integrate updated origin/main with local changes first.\n' >&2; exit 1; }
    git merge --ff-only origin/main
fi
venv="${PINCER_VENV:-$root/.venv}"
if [[ -x "$venv/bin/python" ]]; then
    runner="$venv/bin/python"
elif [[ -f "$venv/Scripts/python.exe" ]]; then
    runner="$venv/Scripts/python.exe"
else
    printf 'ERROR: Phase 3 requires the configured scientific environment (PINCER_VENV).\n' >&2
    exit 1
fi
export PINCER_INSTALL_EXTRAS="${PINCER_INSTALL_EXTRAS:-test,campaign,advanced,documents}"
export PINCER_COMMIT_MESSAGE="${PINCER_COMMIT_MESSAGE:-feat(phase3): add audited spin surfaces, microsolvation, analytic kinetics and native EGNN}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
printf '[phase3] Running all 40 numerical grids and validating all completed native stage receipts.\n'
"$runner" scripts/run_advanced_campaign.py --stage all
printf '[phase3] Entering full tests, format checks, conventional commit and verified push.\n'
exec bash scripts/deploy.sh
