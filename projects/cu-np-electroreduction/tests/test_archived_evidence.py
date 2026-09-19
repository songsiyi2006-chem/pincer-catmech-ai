"""Checks of ACTUAL archived pilot evidence; does not execute quantum chemistry.

The molecular pilot is neither a Cu catalyst calculation nor experimental data.
"""
import hashlib
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from science import HARTREE_EV, parse_xtb


class ArchivedEvidenceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = {}
        for path in sorted((ROOT / "data/raw").glob("pilot_00[12]/*/record.json")):
            record = json.loads(path.read_text(encoding="utf-8"))
            if record["calculation_id"] in cls.records:
                raise AssertionError("Duplicate archived calculation identifier")
            cls.records[record["calculation_id"]] = (path.parent, record)
        cls.summary = json.loads((ROOT / "results/pilot_summary.json").read_text(encoding="utf-8"))
        cls.config = json.loads((ROOT / "config/pilot.json").read_text(encoding="utf-8"))

    def test_all_twenty_native_records_are_reparseable_and_unmodified(self):
        self.assertEqual(len(self.records), 20)
        self.assertEqual(self.summary["native_calls"], 20)
        for identifier, (folder, record) in self.records.items():
            with self.subTest(calculation=identifier):
                self.assertEqual(record["status"], "completed")
                self.assertEqual(record["exit_code"], 0)
                for name, digest in record["output_hashes"].items():
                    self.assertEqual(hashlib.sha256((folder / name).read_bytes()).hexdigest(), digest)
                raw = "\n".join((folder / name).read_text(encoding="utf-8", errors="replace") for name in ("stdout.log", "stderr.log"))
                parsed = parse_xtb(raw, record["settings"]["task"])
                self.assertAlmostEqual(parsed["energy_hartree"], record["energy_hartree"], places=11)

    def test_charge_pairs_use_identical_nuclei_and_correct_electron_parity(self):
        for row in self.summary["fixed_geometry_charging_energies"]:
            negative_folder, negative = self.records[row["calculation_ids"][0]]
            neutral_folder, neutral = self.records[row["calculation_ids"][1]]
            self.assertEqual((negative_folder / "input.xyz").read_bytes(), (neutral_folder / "input.xyz").read_bytes())
            self.assertEqual(negative["settings"]["charge"], -1)
            self.assertEqual(neutral["settings"]["charge"], 0)
            for folder, record in ((negative_folder, negative), (neutral_folder, neutral)):
                lines = (folder / "input.xyz").read_text().splitlines()
                elements = [line.split()[0] for line in lines[2:]]
                self.assertEqual(len(elements), int(lines[0]))
                counts = {element: elements.count(element) for element in set(elements)}
                self.assertEqual(counts, {"C": 9, "H": 8, "N": 2, "O": 1})
                electrons = 9 * 6 + 8 + 2 * 7 + 8 - record["settings"]["charge"]
                self.assertEqual((electrons - record["settings"]["unpaired"]) % 2, 0)
            independent_difference = (negative["energy_hartree"] - neutral["energy_hartree"]) * HARTREE_EV
            self.assertAlmostEqual(row["delta_E_qminus1_minus_q0_eV"], independent_difference, places=10)

    def test_rejected_geometry_is_retained_and_repair_has_54_positive_vibrations(self):
        limit = self.config["predeclared_tolerances"]["neutral_minimum_imaginary_cm1_threshold"]
        for identifier, accepted in (("P009-neutral-gfn2-gas-hess", False), ("R001-rotor-hess", True)):
            folder, record = self.records[identifier]
            text = (folder / "vibspectrum").read_text(encoding="utf-8")
            vibrations = [float(x) for x in re.findall(r"^\s*\d+\s+[ai]\s+(-?\d+\.\d+)", text, re.M)]
            self.assertEqual(len(vibrations), 3 * 20 - 6)
            self.assertEqual(min(vibrations) >= limit, accepted)
            self.assertEqual(record["minimum_check"] == "passed", accepted)
        repaired = self.records["R001-rotor-hess"][1]["frequencies_cm1"]
        self.assertAlmostEqual(min(repaired), 78.64, places=2)

    def test_method_disagreement_triggers_predeclared_interpretation_stop(self):
        limit = self.config["predeclared_tolerances"]["charging_energy_method_spread_flag_eV"]
        self.assertEqual(self.summary["predeclared_spread_threshold_eV"], limit)
        values = {(x["gfn"], x["solvent"]): x["delta_E_qminus1_minus_q0_eV"] for x in self.summary["fixed_geometry_charging_energies"]}
        for solvent in ("gas", "dmf"):
            spread = abs(values[(1, solvent)] - values[(2, solvent)])
            self.assertAlmostEqual(self.summary["hamiltonian_spread_eV"][solvent], spread, places=10)
            self.assertGreater(spread, limit)
        self.assertTrue(self.summary["method_sensitivity_stop"])

    def test_holdout_is_reserved_and_not_used_by_this_molecular_pilot(self):
        holdout = json.loads((ROOT / "config/prospective_holdout.json").read_text(encoding="utf-8"))
        self.assertEqual(len({x["id"] for x in holdout["holdouts"]}), 3)
        for _, record in self.records.values():
            self.assertEqual(record["scope"], "molecular_method_preflight_only")
            self.assertNotIn(record["calculation_id"], {x["id"] for x in holdout["holdouts"]})
        # This does not validate future scaffold grouping or prospective accuracy.


if __name__ == "__main__":
    unittest.main()
