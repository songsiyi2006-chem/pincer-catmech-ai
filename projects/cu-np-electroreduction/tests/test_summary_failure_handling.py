"""SYNTHETIC FAILURE INJECTIONS into temporary copies of the native archive.

No production logs are changed; synthetic failures are not research observations.
"""
import contextlib
import hashlib
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import summarize_pilot


class SummaryFailureHandling(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cu_pilot_review_")
        self.fixture = Path(self.temp.name)
        for dirname in ("config", "data/raw"):
            shutil.copytree(ROOT / dirname, self.fixture / dirname)

    def tearDown(self):
        self.temp.cleanup()

    def regenerate(self):
        with patch.object(summarize_pilot, "ROOT", self.fixture), contextlib.redirect_stdout(io.StringIO()):
            summarize_pilot.main()
        return json.loads((self.fixture / "results/pilot_summary.json").read_text(encoding="utf-8"))

    def test_failed_native_call_remains_indexed_and_never_supplies_an_energy(self):
        folder = self.fixture / "data/raw/pilot_003/P-SYNTHETIC-FAILURE"
        folder.mkdir(parents=True)
        (folder / "stdout.log").write_text("SYNTHETIC SOFTWARE TEST\nTOTAL ENERGY -999 Eh\n", encoding="utf-8")
        (folder / "stderr.log").write_text("SYNTHETIC SOFTWARE TEST\nabnormal termination of xtb\n", encoding="utf-8")
        record = {
            "calculation_id": "P-SYNTHETIC-FAILURE", "status": "failed", "error": "synthetic nonconvergence",
            "settings": {"task": "sp", "gfn": 2, "solvent": "dmf", "charge": -1, "unpaired": 1},
            "elapsed_seconds": 0.0, "identity": {"input_sha256": "synthetic-fixture-no-input"},
            "output_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir()},
        }
        (folder / "record.json").write_text(json.dumps(record), encoding="utf-8")
        summary = self.regenerate()
        failed = next(x for x in summary["rows"] if x["id"] == "P-SYNTHETIC-FAILURE")
        self.assertEqual(summary["native_calls"], 21)
        self.assertEqual(summary["completed_native_calls"], 20)
        self.assertEqual(failed["status"], "failed")
        self.assertIsNone(failed["energy_hartree"])
        self.assertEqual(len(summary["fixed_geometry_charging_energies"]), 4)
        index = json.loads((self.fixture / "data/raw/index.json").read_text(encoding="utf-8"))
        self.assertIn("P-SYNTHETIC-FAILURE", {x["id"] for x in index})

    def test_sensitivity_gate_reads_configuration(self):
        path = self.fixture / "config/pilot.json"
        config = json.loads(path.read_text(encoding="utf-8"))
        config["predeclared_tolerances"]["charging_energy_method_spread_flag_eV"] = 0.9
        path.write_text(json.dumps(config), encoding="utf-8")
        summary = self.regenerate()
        self.assertEqual(summary["predeclared_spread_threshold_eV"], 0.9)
        self.assertFalse(summary["method_sensitivity_stop"])

    def test_failed_hessian_cannot_validate_downstream_pairs(self):
        path = self.fixture / "data/raw/pilot_002/R001-rotor-hess/record.json"
        record = json.loads(path.read_text(encoding="utf-8"))
        record["status"] = "failed"
        record["error"] = "SYNTHETIC SOFTWARE TEST: native exit failed after spectrum appeared"
        path.write_text(json.dumps(record), encoding="utf-8")
        summary = self.regenerate()
        self.assertEqual(summary["fixed_geometry_charging_energies"], [])

    def test_failure_before_log_creation_remains_indexed(self):
        folder = self.fixture / "data/raw/pilot_003/P-SYNTHETIC-NO-LOG"
        folder.mkdir(parents=True)
        record = {
            "calculation_id": "P-SYNTHETIC-NO-LOG", "status": "failed",
            "error": "SYNTHETIC SOFTWARE TEST: process could not start",
            "settings": {"task": "sp", "gfn": 2, "solvent": "gas", "charge": 0, "unpaired": 0},
            "elapsed_seconds": 0.0, "identity": {"input_sha256": "synthetic-fixture-no-input"},
            "output_hashes": {},
        }
        (folder / "record.json").write_text(json.dumps(record), encoding="utf-8")
        summary = self.regenerate()
        failed = next(x for x in summary["rows"] if x["id"] == "P-SYNTHETIC-NO-LOG")
        self.assertEqual(summary["native_calls"], 21)
        self.assertEqual(failed["status"], "failed")
        self.assertIsNone(failed["energy_hartree"])


if __name__ == "__main__":
    unittest.main()
