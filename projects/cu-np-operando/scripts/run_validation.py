"""Run software/native-record tests, never quantum chemistry, into a new folder."""
import argparse
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import platform
import unittest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise SystemExit('Refusing to overwrite a previous validation record')
    args.output.mkdir(parents=True)
    root = Path(__file__).resolve().parents[1]
    suite = unittest.TestLoader().discover(str(root / 'tests'))
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (args.output / 'tests.txt').write_text(stream.getvalue(), encoding='utf-8')
    summary = {'utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
               'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
               'skipped': len(result.skipped), 'passed': result.wasSuccessful(),
               'scope': 'software and native-record integrity; not experimental validation'}
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    raise SystemExit(not result.wasSuccessful())


if __name__ == '__main__':
    main()
