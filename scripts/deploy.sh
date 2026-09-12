#!/usr/bin/env bash
set -Eeuo pipefail

trap 'code=$?; printf "ERROR: deployment failed at line %s (exit %s).\n" "$LINENO" "$code" >&2; exit "$code"' ERR
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '[deploy] %s\n' "$*"; }

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd -- "$root"
command -v git >/dev/null || fail 'git is required.'
[[ "$(git rev-parse --show-toplevel)" -ef "$root" ]] || fail 'Run from a repository root.'
[[ "$(git symbolic-ref --quiet --short HEAD)" == main ]] || fail 'Checkout main before deployment.'
git remote get-url origin >/dev/null || fail 'Configure the origin remote first.'
[[ ! -e "$(git rev-parse --git-path MERGE_HEAD)" ]] || fail 'Finish the existing merge first.'
[[ ! -d "$(git rev-parse --git-path rebase-merge)" && ! -d "$(git rev-parse --git-path rebase-apply)" ]] || fail 'Finish the existing rebase first.'

log 'Fetching origin before verification.'
git fetch --prune origin
if git show-ref --verify --quiet refs/remotes/origin/main; then
    if ! git rev-parse --verify HEAD >/dev/null 2>&1; then
        fail 'Remote main now exists but local main is unborn; integrate it first.'
    fi
    if ! git merge-base --is-ancestor origin/main HEAD; then
        [[ -z "$(git status --porcelain)" ]] || fail 'Remote main advanced; integrate it with local changes first.'
        git merge --ff-only origin/main
    fi
fi

venv="${PINCER_VENV:-$root/.venv}"
bootstrap="${PYTHON:-python}"
if [[ ! -d "$venv" ]]; then
    log "Creating virtual environment: $venv"
    "$bootstrap" -m venv "$venv"
fi
if [[ -x "$venv/bin/python" ]]; then
    runner="$venv/bin/python"
elif [[ -f "$venv/Scripts/python.exe" ]]; then
    runner="$venv/Scripts/python.exe"
else
    fail "Existing environment has no Python interpreter: $venv"
fi

log 'Installing package and test dependencies.'
"$runner" -m pip install -e '.[test]'
log 'Running the full verification gate.'
if ! "$runner" -m pytest tests/ -v; then
    fail 'Pytest failed; no staging, commit, or push was attempted.'
fi
git diff --check
git diff --cached --check

paths=(src tests docs scripts examples .github README.md pyproject.toml .gitignore .gitattributes)
for path in "${paths[@]}"; do
    if [[ -e "$path" ]] || git ls-files --error-unmatch -- "$path" >/dev/null 2>&1; then
        git add --all -- "$path"
    fi
done
while IFS= read -r -d '' staged; do
    case "$staged" in
        src/*|tests/*|docs/*|scripts/*|examples/*|.github/*|README.md|pyproject.toml|.gitignore|.gitattributes) ;;
        *) fail "Staged path is outside the deployment file set: $staged" ;;
    esac
done < <(git diff --cached --name-only -z)
git diff --cached --check
if ! git diff --cached --quiet; then
    git commit -m 'feat(core): implement quasi-rrho thermodynamics, steric profiling, and bilingual docs'
else
    log 'No staged changes; reusing the verified commit.'
fi
git rev-parse --verify HEAD >/dev/null || fail 'Nothing to push: repository has no commit.'
[[ -z "$(git status --porcelain)" ]] || fail 'Uncommitted paths remain outside the deployment file set.'
log 'Pushing verified main without force.'
git push origin main
[[ "$(git rev-parse HEAD)" == "$(git ls-remote origin refs/heads/main | cut -f1)" ]] || fail 'Remote main differs from local HEAD after push.'
log 'Deployment complete; remote main matches local HEAD.'
