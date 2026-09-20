"""Process-control failure tests: no real child or quantum engine is started."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import psutil

import bounded_runner


class BoundedRunnerChecks(unittest.TestCase):
    def test_disappearing_child_does_not_skip_parent_kill(self):
        child = Mock()
        child.kill.side_effect = psutil.NoSuchProcess(12345)
        parent = Mock()
        parent.children.return_value = [child]
        proc = Mock()
        proc.poll.return_value = 0
        bounded_runner.kill_process_tree(parent, proc)
        parent.kill.assert_called_once()

    def test_launch_exception_is_persisted_as_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "new-run"
            geometry = Path(__file__).resolve().parent / "runs/DFT001-pbe0-neutral/input.xyz"
            argv = ["bounded_runner.py", "--psi4-python", "missing-python", "--geometry", str(geometry), "--method", "pbe0", "--charge", "0", "--output", str(output)]
            with patch("sys.argv", argv), patch.object(bounded_runner.psutil, "virtual_memory", return_value=Mock(available=4 * 1024**3)), patch.object(bounded_runner.subprocess, "Popen", side_effect=OSError("test launch failure")), contextlib.redirect_stdout(io.StringIO()):
                return_code = bounded_runner.main()
            record = json.loads((output / "run_record.json").read_text())
            self.assertEqual(record["status"], "launch_failed")
            self.assertNotEqual(return_code, 0)
            self.assertIsNone(record["exit_code"])
            self.assertIn("test launch failure", record["exception"])


if __name__ == "__main__":
    unittest.main()
