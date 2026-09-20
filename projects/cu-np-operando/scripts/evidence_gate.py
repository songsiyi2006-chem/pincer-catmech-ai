"""Fail-closed research readiness checks; attestations are not scientific proof."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def bound_file(root: Path, relative: str) -> Path:
    path = (root / relative).resolve()
    path.relative_to(root.resolve())
    if not path.is_file():
        raise ValueError(f'Missing artifact: {relative}')
    return path


def check_artifact(root: Path, record: dict) -> list[str]:
    errors = []
    try:
        path = bound_file(root, record['path'])
        if digest(path) != record['sha256']:
            errors.append(f'Hash mismatch: {record["path"]}')
    except (ValueError, KeyError, TypeError, OSError) as exc:
        errors.append(str(exc))
    return errors


REQUIREMENTS = {
    'reaction_identity': 'Verified reactant, partner, product and atom/electron mapping',
    'source_site_coordinates': 'Source-matched periodic Cu structures and cell',
    'potential_calibration': 'Solvent-compatible potential reference and calibration',
    'constant_potential_paths': 'Executed accepted constant-potential competing paths',
    'site_density': 'Measured site density with uncertainty for current conversion',
    'kinetic_parameters': 'Traceable accepted rate parameters, not fixture labels',
    'transport_calibration': 'Geometry, transport coefficients and validation',
    'stability_data': 'Time-resolved stability or declared lifetime uncertainty',
    'independent_validation': 'Prospectively registered independent validation',
}


def evaluate(root: Path, evidence: dict) -> dict:
    """No completed scientific ranking unless every required record is substantiated.

    Machine gates validate completeness/integrity only. A signed review remains
    a human or external scientific attestation, not a consequence of a hash.
    """
    results = {}
    requirements = evidence.get('requirements') if isinstance(evidence, dict) else None
    if not isinstance(requirements, dict):
        requirements = {}
    for key, explanation in REQUIREMENTS.items():
        entry = requirements.get(key, {})
        if not isinstance(entry, dict):
            entry = {}
        errors = []
        if entry.get('status') != 'accepted':
            errors.append('not accepted')
        if entry.get('evidence_kind') not in ('LITERATURE', 'CALCULATION', 'EXPERIMENT'):
            errors.append('missing or nonphysical evidence kind')
        if not all(isinstance(entry.get(k), str) and entry[k].strip()
                   for k in ('scope_review', 'reviewer')):
            errors.append('scope/reviewer attestation absent')
        artifacts = entry.get('artifacts', [])
        if not isinstance(artifacts, list) or not artifacts:
            errors.append('no local supporting artifact')
            artifacts = []
        for item in artifacts:
            errors.extend(check_artifact(root, item))
        results[key] = {'passed': not errors, 'definition': explanation, 'errors': errors}
    return {
        'evidence_manifest_complete': all(v['passed'] for v in results.values()),
        'physical_ranking_ready': False,
        'blocking_reason': 'No validated target-reaction model is implemented; metadata cannot unlock physical ranking',
        'scope': 'completeness gate, not an independent scientific certificate',
        'requirements': results,
    }


def audit_group_split(rows: list[dict]) -> dict:
    """The entire geometry/state family must remain in a single split."""
    family_splits = {}
    errors = []
    for row in rows:
        family, split = row.get('family_id'), row.get('split')
        if not family or split not in ('train', 'validation', 'test'):
            errors.append('row has missing family_id or invalid split')
            continue
        family_splits.setdefault(family, set()).add(split)
    conflicts = {k: sorted(v) for k, v in family_splits.items() if len(v) > 1}
    return {'passed': not errors and not conflicts, 'conflicts': conflicts, 'errors': errors,
            'scope': 'declared family consistency only; chemical equivalence needs separate review'}


def manifest(root: Path) -> dict:
    excluded = {'__pycache__', '.pytest_cache', '.git'}
    records = []
    for path in sorted(root.rglob('*')):
        rel = path.relative_to(root)
        if not path.is_file() or excluded.intersection(rel.parts):
            continue
        if rel.as_posix() == 'results/file_manifest.json' or path.suffix in ('.pyc', '.pyo'):
            continue
        records.append({'path': rel.as_posix(), 'bytes': path.stat().st_size, 'sha256': digest(path)})
    return {'algorithm': 'SHA256', 'scope': 'this project only; excludes caches and self', 'files': records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['gate', 'manifest', 'verify'])
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    results = root / 'results'
    results.mkdir(exist_ok=True)
    if args.action == 'gate':
        result = evaluate(root, json.loads((root / 'config/evidence_requirements.json').read_text(encoding='utf-8')))
        (results / 'production_readiness.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        print(json.dumps(result, indent=2))
        raise SystemExit(0 if result['physical_ranking_ready'] else 2)
    elif args.action == 'manifest':
        result = manifest(root)
        (results / 'file_manifest.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        print(f'Manifested {len(result["files"])} files')
    else:
        expected = json.loads((results / 'file_manifest.json').read_text(encoding='utf-8'))
        errors = [e for record in expected['files'] for e in check_artifact(root, record)]
        actual_paths = {r['path'] for r in manifest(root)['files']}
        expected_paths = {r['path'] for r in expected['files']}
        errors += [f'Unexpected unmanifested file: {p}' for p in sorted(actual_paths - expected_paths)]
        print(json.dumps({'passed': not errors, 'file_count': len(expected_paths), 'errors': errors}, indent=2))
        raise SystemExit(bool(errors))


if __name__ == '__main__':
    main()
