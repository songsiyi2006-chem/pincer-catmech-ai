"""Archive validation only: no Psi4 execution and no invented science labels."""
import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from summarize_dft import HARTREE_TO_EV, sha256, summarize, validate_pair_protocol

ROOT = Path(__file__).resolve().parent


class DFTArchiveChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copytree(ROOT / "runs", self.root / "runs", ignore=shutil.ignore_patterns("scratch", "psi.*"))

    def summary(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return summarize(self.root)

    def refresh_result_digest(self, folder):
        record_path = folder / "run_record.json"
        record = json.loads(record_path.read_text())
        record["files"]["result.json"] = sha256(folder / "result.json")
        record_path.write_text(json.dumps(record))

    def completed(self):
        return [p.parent for p in (self.root / "runs").glob("*/result.json") if json.loads(p.read_text())["status"] == "completed"]

    def test_native_hash_energy_electron_and_protocol_validation(self):
        result = self.summary()
        self.assertGreater(result["completed_quantum_runs"], 0)

    def test_actual_pair_numeric_sign_and_units(self):
        result = self.summary()
        if not result["energy_pairs"]:
            self.skipTest("No completed neutral/anion pair; physical difference remains unavailable")
        for pair in result["energy_pairs"]:
            neutral = next(r for r in result["runs"] if r["run_id"] == pair["neutral_run_id"])
            anion = next(r for r in result["runs"] if r["run_id"] == pair["anion_run_id"])
            independently_computed = (neutral["energy_hartree"] - anion["energy_hartree"]) * HARTREE_TO_EV
            self.assertAlmostEqual(pair["electronic_attachment_energy_neutral_minus_anion_eV"], independently_computed, places=12)
            self.assertAlmostEqual(pair["delta_E_anion_minus_neutral_eV"], -independently_computed, places=12)

    def test_tampered_native_output_is_rejected(self):
        folder = self.completed()[0]
        with (folder / "psi4.out").open("a") as f:
            f.write("\nCORRUPTED TEST COPY\n")
        with self.assertRaisesRegex(ValueError, "Modified archived output"):
            self.summary()

    def test_electron_count_mutation_cannot_supply_a_label(self):
        folder = self.completed()[0]
        result_path = folder / "result.json"
        result = json.loads(result_path.read_text())
        result["nalpha"] += 1
        result_path.write_text(json.dumps(result))
        self.refresh_result_digest(folder)
        with self.assertRaisesRegex(ValueError, "Electron count mismatch"):
            self.summary()

    def test_spin_inconsistency_is_rejected(self):
        folder = self.completed()[0]
        result_path = folder / "result.json"
        result = json.loads(result_path.read_text())
        result["s2_ideal"] = 99.0
        result_path.write_text(json.dumps(result))
        self.refresh_result_digest(folder)
        with self.assertRaisesRegex(ValueError, "Ideal spin expectation"):
            self.summary()

    def test_nan_energy_cannot_supply_a_label(self):
        folder = self.completed()[0]
        result_path = folder / "result.json"
        result = json.loads(result_path.read_text())
        result["energy_hartree"] = float("nan")
        result_path.write_text(json.dumps(result))
        self.refresh_result_digest(folder)
        with self.assertRaisesRegex(ValueError, "Nonfinite result field"):
            self.summary()

    def test_settings_conflict_is_rejected(self):
        folder = self.completed()[0]
        result_path = folder / "result.json"
        result = json.loads(result_path.read_text())
        result["settings"]["basis"] = "unrecorded-basis"
        result_path.write_text(json.dumps(result))
        self.refresh_result_digest(folder)
        with self.assertRaisesRegex(ValueError, "protocol fields"):
            self.summary()

    def test_embedded_script_hash_conflict_is_rejected(self):
        folder = self.completed()[0]
        result_path = folder / "result.json"
        result = json.loads(result_path.read_text())
        result["child_script_sha256"] = "wrong"
        result_path.write_text(json.dumps(result))
        self.refresh_result_digest(folder)
        with self.assertRaisesRegex(ValueError, "child-script hash"):
            self.summary()

    def test_missing_external_preflight_is_not_invented(self):
        self.assertFalse((self.root / "memory_preflight_rejection.json").exists())
        result = self.summary()
        actual = sum(r["status"] == "resource_preflight_not_launched" for r in result["runs"])
        self.assertEqual(result["preflight_rejections_not_counted_as_quantum_runs"], actual)

    def test_duplicate_success_is_rejected(self):
        shutil.copytree(self.completed()[0], self.root / "runs" / "ZZZ-ambiguous-success-copy")
        with self.assertRaisesRegex(ValueError, "Ambiguous duplicate"):
            self.summary()

    def test_failed_run_cannot_supply_a_charge_pair(self):
        folder = self.completed()[0]
        path = folder / "run_record.json"
        record = json.loads(path.read_text())
        record["status"] = "timeout"
        path.write_text(json.dumps(record))
        result = self.summary()
        self.assertGreaterEqual(result["failed_or_terminated_runs"], 1)
        self.assertTrue(all(folder.name not in [p["neutral_run_id"], p["anion_run_id"]] for p in result["energy_pairs"]))

    def test_pair_geometry_basis_and_grid_mismatch_rejected(self):
        actual = json.loads((self.completed()[0] / "result.json").read_text())
        # Copies are validation fixtures only; never used to create energy labels.
        for mutation, message in [
            (lambda r: r.update(input_sha256="different"), "geometry hash"),
            (lambda r: r["settings"].update(basis="different"), "energy-pair setting"),
            (lambda r: r["psi4_options"].update(dft_radial_points=999), "SCF/DF/grid protocol"),
            (lambda r: r.update(basis_source_sha256="different"), "basis source content"),
            (lambda r: r.update(psi4_version="different"), "Psi4 versions"),
        ]:
            with self.subTest(message=message):
                changed = copy.deepcopy(actual)
                mutation(changed)
                with self.assertRaisesRegex(ValueError, message):
                    validate_pair_protocol(actual, changed)


if __name__ == "__main__":
    unittest.main()
