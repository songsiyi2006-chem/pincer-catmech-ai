"""Expose native-record checks without importing or running Psi4."""
from pathlib import Path
import sys

ARCHIVE = Path(__file__).resolve().parents[1] / 'data' / 'dft_pilot'
sys.path.insert(0, str(ARCHIVE))
from test_dft_archive import DFTArchiveChecks  # noqa: E402,F401
from test_bounded_runner import BoundedRunnerChecks  # noqa: E402,F401
