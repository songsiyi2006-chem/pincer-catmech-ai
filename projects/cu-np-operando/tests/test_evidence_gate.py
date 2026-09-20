import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from evidence_gate import REQUIREMENTS, audit_group_split, check_artifact, evaluate


class EvidenceGateTests(unittest.TestCase):
    def test_no_evidence_never_ready(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertFalse(evaluate(Path(folder), {})['physical_ranking_ready'])

    def test_synthetic_cannot_unlock_physical_prediction(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'fixture.json').write_text('{}')
            record = {'status': 'accepted', 'evidence_kind': 'SYNTHETIC', 'reviewer': 'test',
                      'scope_review': 'test only', 'artifacts': [
                          {'path': 'fixture.json', 'sha256': hashlib.sha256(b'{}').hexdigest()}]}
            evidence = {'requirements': {k: record for k in REQUIREMENTS}}
            self.assertFalse(evaluate(root, evidence)['physical_ranking_ready'])

    def test_artifact_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'data').write_bytes(b'original')
            record = {'path': 'data', 'sha256': hashlib.sha256(b'original').hexdigest()}
            self.assertEqual(check_artifact(root, record), [])
            (root / 'data').write_bytes(b'changed')
            self.assertTrue(check_artifact(root, record))

    def test_escape_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertTrue(check_artifact(Path(folder), {'path': '../outside', 'sha256': 'x'}))

    def test_family_leakage_detected(self):
        rows = [{'family_id': 'substrate-a', 'split': 'train', 'potential': -1},
                {'family_id': 'substrate-a', 'split': 'test', 'potential': -2}]
        self.assertFalse(audit_group_split(rows)['passed'])

    def test_grouped_split(self):
        rows = [{'family_id': 'a', 'split': 'train'}, {'family_id': 'b', 'split': 'test'}]
        self.assertTrue(audit_group_split(rows)['passed'])

    def test_relabelled_fixture_cannot_unlock_target_model(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            content = b'SYNTHETIC_UNIT_FIXTURE'
            (root / 'fixture').write_bytes(content)
            entry = {'status': 'accepted', 'evidence_kind': 'CALCULATION',
                     'reviewer': 'claimed', 'scope_review': 'claimed',
                     'artifacts': [{'path': 'fixture', 'sha256': hashlib.sha256(content).hexdigest()}]}
            result = evaluate(root, {'requirements': {k: entry for k in REQUIREMENTS}})
            self.assertTrue(result['evidence_manifest_complete'])
            self.assertFalse(result['physical_ranking_ready'])

    def test_malformed_requirements_fail_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            for bad in (None, {'requirements': None}, {'requirements': []},
                        {'requirements': {'reaction_identity': None}}):
                self.assertFalse(evaluate(Path(folder), bad)['physical_ranking_ready'])


if __name__ == '__main__':
    unittest.main()
