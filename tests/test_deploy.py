"""Exercise deployment side effects using disposable local Git remotes.

The Python shim controls the shell gate outcome; actual scientific tests run
under pytest separately. No test contacts GitHub or modifies the real checkout.
"""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


def _git(cwd, *args):
    return subprocess.run(
        [shutil.which("git"), "-C", str(cwd), *args], check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.fixture
def deployment(request):
    # Own the directory lifecycle. Repeated pytest TempPathFactory operations
    # can break Windows/MSYS process creation even for `bash -c true`.
    temporary = tempfile.TemporaryDirectory(prefix="pincer-deploy-")
    request.addfinalizer(temporary.cleanup)
    tmp_path = Path(temporary.name)
    bash = shutil.which("bash")
    windows_bash = Path("C:/Program Files/Git/bin/bash.exe")
    if os.name == "nt" and windows_bash.exists():
        bash = str(windows_bash)
    if not bash or not shutil.which("git"):
        pytest.skip("Bash and git are needed for deployment integration tests")
    remote = tmp_path / "remote.git"
    checkout = tmp_path / "checkout with spaces"
    remote.mkdir()
    checkout.mkdir()
    _git(remote, "init", "--bare", "--initial-branch=main")
    _git(checkout, "init", "--initial-branch=main")
    _git(checkout, "config", "user.name", "Deployment Test")
    _git(checkout, "config", "user.email", "deployment@example.invalid")
    _git(checkout, "config", "commit.gpgsign", "false")
    _git(checkout, "config", "core.autocrlf", "false")
    _git(checkout, "remote", "add", "origin", remote.as_posix())
    scripts = checkout / "scripts"
    scripts.mkdir()
    source = Path(__file__).resolve().parents[1] / "scripts" / "deploy.sh"
    shutil.copyfile(source, scripts / "deploy.sh")
    (checkout / "README.md").write_text("# Integration fixture\n", encoding="utf-8", newline="\n")
    venv = tmp_path / "controlled venv"
    runner = venv / "bin" / "python"
    runner.parent.mkdir(parents=True)
    runner.write_text(
        '#!/usr/bin/env bash\n'
        'set -eu\n'
        'if [[ "$*" == "-m pip install -e .[test]" ]]; then exit 0; fi\n'
        'if [[ "$*" == "-m pytest tests/ -v" ]]; then\n'
        '  printf "CONTROLLED_PYTEST_GATE=%s\\n" "${GATE_EXIT:-0}"\n'
        '  exit "${GATE_EXIT:-0}"\n'
        'fi\n'
        'exit 91\n', encoding="utf-8", newline="\n"
    )
    runner.chmod(0o755)
    env = dict(os.environ, PINCER_VENV=venv.as_posix(), GIT_TERMINAL_PROMPT="0")

    def run(exit_code=0):
        command = [bash, "--noprofile", "--norc", (scripts / "deploy.sh").as_posix()]
        return subprocess.run(
            command,
            env=dict(env, GATE_EXIT=str(exit_code)), capture_output=True,
            text=True, encoding="utf-8", errors="replace", timeout=60,
        )

    return checkout, remote, run


def test_failed_tests_never_stage_commit_or_push(deployment):
    checkout, remote, run = deployment
    result = run(7)
    assert result.returncode != 0
    assert "Pytest failed; no staging, commit, or push" in result.stderr
    assert not _git(checkout, "ls-files")
    assert not _git(remote, "for-each-ref", "refs/heads/main")


def test_successful_gate_pushes_and_second_run_is_idempotent(deployment):
    checkout, remote, run = deployment
    first = run()
    assert first.returncode == 0, first.stdout + first.stderr
    head = _git(checkout, "rev-parse", "HEAD")
    assert head == _git(remote, "rev-parse", "main")
    assert "CONTROLLED_PYTEST_GATE=0" in first.stdout
    assert _git(checkout, "log", "-1", "--format=%s") == (
        "feat(core): implement quasi-rrho thermodynamics, steric profiling, and bilingual docs"
    )
    second = run()
    assert second.returncode == 0, second.stdout + second.stderr
    assert _git(checkout, "rev-parse", "HEAD") == head


def test_staged_outside_project_is_not_committed(deployment):
    checkout, remote, run = deployment
    (checkout / "private.txt").write_text("unrelated\n", encoding="utf-8", newline="\n")
    _git(checkout, "add", "private.txt")
    result = run()
    assert result.returncode != 0
    assert "outside the deployment file set" in result.stderr
    assert not _git(remote, "for-each-ref", "refs/heads/main")


def test_wrong_branch_stops_before_test_gate(deployment):
    checkout, remote, run = deployment
    _git(checkout, "symbolic-ref", "HEAD", "refs/heads/topic")
    result = run()
    assert result.returncode != 0
    assert "Checkout main" in result.stderr
    assert "CONTROLLED_PYTEST" not in result.stdout
    assert not _git(remote, "for-each-ref", "refs/heads/main")
