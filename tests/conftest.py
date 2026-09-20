"""Portable temporary paths for scientific and shell integration fixtures."""

import os
from pathlib import Path
import tempfile

import pytest


if os.name == "nt":
    @pytest.fixture
    def tmp_path():
        """Keep Windows process startup independent of TempPathFactory state.

        On the tested Windows/Python host, repeated pytest-managed paths cause
        later CreateProcess calls to fail with WinError 5, even for bash -c true.
        Standard-library temporary directories preserve per-test isolation and
        deterministic cleanup without changing subprocess errors or test results.
        """
        with tempfile.TemporaryDirectory(prefix="pincer-test-") as directory:
            # Windows runners may expose TEMP through an 8.3 alias.
            # Match production roots, which are already canonicalized.
            yield Path(directory).resolve()
