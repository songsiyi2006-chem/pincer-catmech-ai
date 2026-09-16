"""Portable, standard-library-only ORCA/Slurm bridge. Not a kinetic certificate.

This file is copied into every bundle, so cluster execution needs Python >=3.10
and site-provided ORCA/MPI, but no installation of this Python package.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
import zipfile


ELEMENTS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe '
            'Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In '
            'Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf '
            'Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am '
            'Cm Bk Cf Es Fm Md No Lr').split()
Z = {symbol: index + 1 for index, symbol in enumerate(ELEMENTS)}
SAFE_ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}\Z')
SHA = re.compile(r'[0-9a-f]{64}\Z')
PROTOCOL = {
    'engine': 'ORCA', 'supported_major_version': 6, 'task': 'single_point',
    'method': 'r2SCAN-3c', 'basis': 'method-internal def2-mTZVPP',
    'solvation_model': 'CPCM', 'solvent': 'toluene', 'scf': 'TightSCF',
    'max_scf_iterations': 500, 'restricted_singlet': True,
    'reaction_reference_temperature_K': 383.15,
    'temperature_used_in_electronic_energy': False,
    'dielectric_temperature_correction': 'not_applied_ORCA_builtin_solvent',
    'purpose': 'geometry_vertical_spin_diagnostic',
    'experimental_identity_verified': False,
    'physical_kinetic_certificate': False,
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest_object(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def load_json(path):
    def pairs(values):
        out = {}
        for key, val in values:
            if key in out:
                raise ValueError(f'duplicate JSON key: {key}')
            out[key] = val
        return out
    def bad_constant(value):
        raise ValueError(f'nonfinite JSON constant: {value}')
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      object_pairs_hook=pairs, parse_constant=bad_constant)


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                    allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def integer(value, name, lower=1, upper=1000000):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be integer in [{lower}, {upper}]')
    return value


def relative_file(root, name):
    """No extraction, absolute path, symlink, or traversal through the bundle."""
    if not isinstance(name, str) or '\\' in name or ':' in name:
        raise ValueError('path must be a portable relative POSIX path')
    rel = PurePosixPath(name)
    if rel.is_absolute() or not rel.parts or any(p in ('', '.', '..') for p in name.split('/')):
        raise ValueError('unsafe relative path')
    root = Path(root).resolve()
    candidate = root
    for part in rel.parts:
        candidate /= part
        if candidate.is_symlink():
            raise ValueError('symlinks are not allowed in evidence')
    if not candidate.resolve().is_relative_to(root):
        raise ValueError('path escapes source root')
    if not candidate.is_file():
        raise ValueError(f'missing source file: {name}')
    return candidate


def read_xyz(path):
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    if len(lines) < 3:
        raise ValueError('truncated XYZ')
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError('invalid XYZ atom count') from exc
    if not 1 <= count <= 10000 or len(lines) != count + 2:
        raise ValueError('XYZ must have exactly one complete frame')
    atoms = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4 or fields[0] not in Z:
            raise ValueError('unsupported element or malformed XYZ row')
        xyz = tuple(float(x) for x in fields[1:])
        if not all(math.isfinite(x) and abs(x) < 1e6 for x in xyz):
            raise ValueError('nonfinite or implausible XYZ coordinate')
        atoms.append((fields[0], *xyz))
    for i, first in enumerate(atoms):
        for second in atoms[i + 1:]:
            if math.dist(first[1:], second[1:]) < .25:
                raise ValueError('overlapping atoms below 0.25 Angstrom')
    return atoms


def validate_electrons(atoms, charge, multiplicity):
    integer(charge, 'charge', -20, 20)
    integer(multiplicity, 'multiplicity', 1, 31)
    electrons = sum(Z[a[0]] for a in atoms) - charge
    if electrons <= 0 or multiplicity - 1 > electrons or (electrons - multiplicity + 1) % 2:
        raise ValueError('electron count / charge / multiplicity parity mismatch')
    return electrons


def validate_metadata(atoms, metadata, charge):
    if type(metadata.get('charge')) is not int or metadata['charge'] != charge:
        raise ValueError('explicit charge must agree with original metadata')
    formula = metadata.get('formula')
    if not isinstance(formula, str) or not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', formula):
        raise ValueError('original metadata must include a plain molecular formula')
    counts = Counter()
    for symbol, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        if symbol not in Z or (count and int(count) < 1):
            raise ValueError('invalid metadata formula')
        counts[symbol] += int(count or '1')
    if counts != Counter(a[0] for a in atoms):
        raise ValueError('metadata formula and XYZ composition mismatch')


def structure_provenance(row, metadata):
    origin = row.get('structure_origin', 'designed_prototype')
    if origin not in ('designed_prototype', 'public_crystal_derived_candidate'):
        raise ValueError('unsupported structure_origin')
    provenance = metadata.get('provenance', {})
    if origin == 'public_crystal_derived_candidate':
        if not isinstance(provenance, dict) or any(
                not isinstance(provenance.get(key), str) or not provenance[key].strip()
                for key in ('source_url', 'source_identifier', 'structure_derivation')):
            raise ValueError('public crystal origin requires explicit source and derivation provenance')
        if not re.match(r'https?://[^/\s]+', provenance['source_url']):
            raise ValueError('public source_url must be an HTTP(S) source citation')
        label = 'public crystal-derived candidate; solution identity and ground spin unverified'
    else:
        label = 'computed design geometry; not experiment-matched'
    return origin, provenance, label


def resources(nprocs, memory_mb, maxcore_mb, walltime, concurrency):
    integer(nprocs, 'nprocs', 1, 256)
    integer(memory_mb, 'memory_mb', 512, 2000000)
    integer(maxcore_mb, 'maxcore_mb', 64, 2000000)
    integer(concurrency, 'concurrency', 1, 128)
    if nprocs * maxcore_mb > .75 * memory_mb:
        raise ValueError('MaxCore * nprocs exceeds 75% of requested node memory')
    if not re.fullmatch(r'\d{2,3}:[0-5]\d:[0-5]\d', walltime) or walltime == '00:00:00':
        raise ValueError('walltime must be nonzero HH:MM:SS')
    return dict(nprocs=nprocs, memory_mb=memory_mb, maxcore_mb=maxcore_mb,
                walltime=walltime, concurrency=concurrency, nodes=1, cpus_per_task=1)


def orca_input(atoms, charge, multiplicity, resource):
    reference = 'RKS' if multiplicity == 1 else 'UKS'
    body = '\n'.join(f'{a[0]:2s} {a[1]:.14f} {a[2]:.14f} {a[3]:.14f}' for a in atoms)
    return (f'! r2SCAN-3c CPCM(toluene) TightSCF SP {reference}\n'
            f'%pal nprocs {resource["nprocs"]} end\n'
            f'%maxcore {resource["maxcore_mb"]}\n'
            '%scf MaxIter 500 end\n'
            f'* xyz {charge} {multiplicity}\n{body}\n*\n')


def prepare(selection_path, source_root, output, resource):
    selection = load_json(selection_path)
    if selection.get('schema') != 'pincer_hpc_selection_v1' or not selection.get('records'):
        raise ValueError('selection schema / nonempty records required')
    resource = resources(**{k: resource[k] for k in
                           ('nprocs', 'memory_mb', 'maxcore_mb', 'walltime', 'concurrency')})
    records, ids = [], set()
    for row in selection['records']:
        label = row.get('id')
        if not isinstance(label, str) or not SAFE_ID.fullmatch(label) or label in ids:
            raise ValueError('unsafe or duplicate task id')
        ids.add(label)
        xyz = relative_file(source_root, row['xyz'])
        metadata = relative_file(source_root, row['metadata'])
        meta = load_json(metadata)
        atoms = read_xyz(xyz)
        charge, multiplicity = row['charge'], row['multiplicity']
        electrons = validate_electrons(atoms, charge, multiplicity)
        validate_metadata(atoms, meta, charge)
        structure_provenance(row, meta)
        records.append((row, xyz, metadata, atoms, electrons))
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError('immutable bundle destination already exists')
    output.parent.mkdir(parents=True, exist_ok=True)
    # Destination is a new sibling, published only after successful verification.
    staging = output.parent / (output.name + '.preparing-' + uuid.uuid4().hex)
    staging.mkdir()
    tasks = []
    for index, (row, xyz, metadata, atoms, electrons) in enumerate(records):
        folder = staging / 'tasks' / row['id']
        folder.mkdir(parents=True)
        shutil.copyfile(xyz, folder / 'geometry.xyz')
        shutil.copyfile(metadata, folder / 'metadata.json')
        inp = orca_input(atoms, row['charge'], row['multiplicity'], resource)
        (folder / 'job.inp').write_text(inp, encoding='utf-8', newline='\n')
        origin, provenance, provenance_label = structure_provenance(row, load_json(metadata))
        tasks.append(dict(index=index, id=row['id'], charge=row['charge'],
                          multiplicity=row['multiplicity'], electron_count=electrons,
                          atom_count=len(atoms), composition=dict(Counter(a[0] for a in atoms)),
                          source_xyz=row['xyz'], source_metadata=row['metadata'],
                          geometry_sha256=sha256(folder / 'geometry.xyz'),
                          metadata_sha256=sha256(folder / 'metadata.json'),
                          input_sha256=sha256(folder / 'job.inp'),
                          experimental_identity_verified=False,
                          solution_geometry_verified=False, ground_spin_verified=False,
                          structure_origin=origin, source_provenance=provenance,
                          provenance=provenance_label))
    shutil.copyfile(__file__, staging / 'bridge.py')
    r = resource
    batch = (f'#!/usr/bin/env bash\n#SBATCH --nodes=1\n#SBATCH --ntasks={r["nprocs"]}\n'
             '#SBATCH --cpus-per-task=1\n'
             f'#SBATCH --mem={r["memory_mb"]}M\n#SBATCH --time={r["walltime"]}\n'
             '#SBATCH --output=slurm-%A_%a.out\n#SBATCH --error=slurm-%A_%a.err\n'
             'set -euo pipefail\n: "${HPC_PYTHON:?set absolute site Python 3.10+ path}"\n'
             'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\n'
             'exec "$HPC_PYTHON" bridge.py run-task --bundle . --index "${SLURM_ARRAY_TASK_ID:?}"\n')
    (staging / 'array.sbatch').write_text(batch, encoding='utf-8', newline='\n')
    readme = ['# Offline ORCA diagnostic bundle', '',
              'Status: prepared_not_submitted. No ORCA calculation is certified by preparation.', '',
              'Protocol: r2SCAN-3c / CPCM(toluene) / TightSCF / single point; ORCA 6.x.',
              '383.15 K is a reaction reference, not a computed thermochemical correction.',
              'No TS, MECP, solution active species, or kinetic certificate is created.', '',
              'Run `python bridge.py verify --bundle .` after transfer.',
              'Site inputs: HPC_ACCOUNT, HPC_PARTITION, ORCA_EXE, HPC_PYTHON, HPC_SCRATCH.',
              'A licensed site installation and administrator-reviewed MPI environment are required.',
              'Use `python bridge.py submit --bundle . --mode pilot` to inspect the command.',
              'Only adding `--execute` calls sbatch; start with one pilot and inspect native output.',
              'Production requires a bound human pilot_review.json; see docs/HPC_BRIDGE_ZH.md in the repository.', '',
              '## Per-structure origin', '']
    for task in tasks:
        readme += [f'- {task["id"]}: {task["structure_origin"]}; {task["provenance"]}.',
                   '  Source metadata (not independently authenticated by this bridge): ' +
                   json.dumps(task['source_provenance'], ensure_ascii=False)]
    (staging / 'README.md').write_text('\n'.join(readme) + '\n', encoding='utf-8', newline='\n')
    files = {p.relative_to(staging).as_posix(): sha256(p)
             for p in sorted(staging.rglob('*')) if p.is_file()}
    manifest = dict(schema='pincer_hpc_bundle_v1', prepared_utc=datetime.now(timezone.utc).isoformat(),
                    status='prepared_not_submitted', protocol=PROTOCOL,
                    protocol_sha256=digest_object(PROTOCOL), resources=resource, tasks=tasks,
                    selection_sha256=sha256(selection_path), files=files,
                    physical_kinetic_certificate=False)
    dump(staging / 'manifest.json', manifest)
    verify_bundle(staging)
    staging.rename(output)
    return manifest


def verify_bundle(bundle):
    bundle = Path(bundle).resolve()
    manifest = load_json(relative_file(bundle, 'manifest.json'))
    if manifest.get('schema') != 'pincer_hpc_bundle_v1' or manifest.get('protocol') != PROTOCOL:
        raise ValueError('unsupported manifest schema / protocol')
    if manifest.get('protocol_sha256') != digest_object(PROTOCOL):
        raise ValueError('protocol hash mismatch')
    r = manifest['resources']
    if r != resources(**{k: r[k] for k in ('nprocs', 'memory_mb', 'maxcore_mb', 'walltime', 'concurrency')}):
        raise ValueError('resource declaration mismatch')
    if not manifest.get('tasks') or not manifest.get('files'):
        raise ValueError('missing tasks or evidence index')
    for name, expected in manifest['files'].items():
        if not SHA.fullmatch(expected) or sha256(relative_file(bundle, name)) != expected:
            raise ValueError(f'bundle hash mismatch: {name}')
    expected_files = {'bridge.py', 'array.sbatch', 'README.md'}
    ids = set()
    for index, task in enumerate(manifest['tasks']):
        if task['index'] != index or not SAFE_ID.fullmatch(task['id']) or task['id'] in ids:
            raise ValueError('task index or identity mismatch')
        ids.add(task['id'])
        prefix = 'tasks/' + task['id'] + '/'
        expected_files.update(prefix + name for name in ('job.inp', 'geometry.xyz', 'metadata.json'))
        atoms = read_xyz(relative_file(bundle, prefix + 'geometry.xyz'))
        ne = validate_electrons(atoms, task['charge'], task['multiplicity'])
        if task['electron_count'] != ne or task['atom_count'] != len(atoms):
            raise ValueError('atom or electron count mismatch')
        if task['composition'] != dict(Counter(a[0] for a in atoms)):
            raise ValueError('composition mismatch')
        metadata = load_json(relative_file(bundle, prefix + 'metadata.json'))
        validate_metadata(atoms, metadata, task['charge'])
        origin, provenance, label = structure_provenance(task, metadata)
        if (task.get('source_provenance') != provenance or task.get('provenance') != label
                or any(task.get(flag) is not False for flag in
                       ('experimental_identity_verified', 'solution_geometry_verified', 'ground_spin_verified'))):
            raise ValueError('structure provenance or evidence boundary mismatch')
        for name, key in [('geometry.xyz', 'geometry_sha256'), ('metadata.json', 'metadata_sha256'),
                          ('job.inp', 'input_sha256')]:
            if manifest['files'].get(prefix + name) != task[key]:
                raise ValueError('task hash binding mismatch')
        expected_input = orca_input(atoms, task['charge'], task['multiplicity'], r)
        if relative_file(bundle, prefix + 'job.inp').read_text(encoding='utf-8') != expected_input:
            raise ValueError('input not generated by declared protocol')
    if set(manifest['files']) != expected_files:
        raise ValueError('unexpected or omitted bundle evidence file')
    return manifest


def parse_orca_output(text, multiplicity):
    """Single-point diagnostics only; no missing field becomes a success."""
    integer(multiplicity, 'multiplicity', 1, 31)
    versions = re.findall(r'^\s*Program Version\s+(\S+)', text, re.MULTILINE)
    energies = re.findall(r'^\s*FINAL SINGLE POINT ENERGY\s+(\S+)\s*$', text, re.MULTILINE)
    s2s = re.findall(r'^\s*Expectation value of <S\*\*2>\s*:\s*(\S+)', text, re.MULTILINE)
    def number(value):
        if value is None:
            return None
        result = float(value.replace('D', 'E').replace('d', 'e'))
        if not math.isfinite(result):
            raise ValueError('nonfinite ORCA numeric field')
        return result
    energy = number(energies[-1]) if energies else None
    s2 = number(s2s[-1]) if s2s else None
    expected = (multiplicity - 1) * (multiplicity + 1) / 4
    normal = len(re.findall(r'^\s*\*+\s*ORCA TERMINATED NORMALLY\s*\*+\s*$', text, re.MULTILINE)) == 1
    scf = bool(re.search(r'^\s*\**\s*SCF CONVERGED AFTER\s+\d+\s+CYCLES', text, re.MULTILINE))
    failure = bool(re.search(r'SCF NOT CONVERGED|SCF DID NOT CONVERGE|ORCA finished by error termination', text))
    version = versions[0] if len(versions) == 1 else None
    correct_version = bool(version and re.fullmatch(r'6\.\d+(?:\.\d+)?(?:[-.][A-Za-z0-9]+)*', version))
    ordered = bool(energies and text.rfind('SCF CONVERGED AFTER') < text.rfind('FINAL SINGLE POINT ENERGY')
                   < text.rfind('ORCA TERMINATED NORMALLY'))
    converged = normal and scf and not failure and len(energies) == 1 and energy is not None and ordered and correct_version
    quality = 'unknown' if s2 is None else ('pass_diagnostic' if abs(s2 - expected) <= .1 else 'reject_contaminated')
    return dict(normal_termination=normal, scf_converged=scf and not failure,
                fixture_marker_detected='SOFTWARE FIXTURE' in text,
                task_converged=converged, engine_version=version,
                engine_version_sha256=hashlib.sha256(version.encode()).hexdigest() if version else None,
                energy_hartree=energy, spin_squared=s2, expected_spin_squared=expected,
                spin_quality=quality, imaginary_frequency_count=None,
                frequency_status='not_calculated_single_point',
                eligible_for_vertical_spin_diagnostic=converged and quality == 'pass_diagnostic',
                certified_minimum=False, certified_transition_state=False,
                certified_mecp=False, physical_kinetic_certificate=False)


def checked_env_path(env, name, directory=False):
    value = env.get(name)
    if not value or not Path(value).is_absolute():
        raise ValueError(f'{name} requires an explicit absolute site path')
    path = Path(value).resolve()
    if not (path.is_dir() if directory else path.is_file()):
        raise ValueError(f'{name} path does not exist')
    return path


def validate_allocation(env, r):
    for name in ('SLURM_JOB_ID', 'SLURM_NTASKS', 'SLURM_JOB_NUM_NODES', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_NODE'):
        if name not in env:
            raise ValueError(f'missing scheduler allocation: {name}')
    if int(env['SLURM_JOB_NUM_NODES']) != 1 or int(env['SLURM_CPUS_PER_TASK']) != 1:
        raise ValueError('bridge supports one-node MPI with one CPU per rank only')
    if int(env['SLURM_NTASKS']) < r['nprocs'] or int(env['SLURM_MEM_PER_NODE']) < r['memory_mb']:
        raise ValueError('scheduler allocation would oversubscribe CPU or memory')


def run_task(bundle, index, env=None):
    env = dict(os.environ if env is None else env)
    bundle = Path(bundle).resolve()
    manifest = verify_bundle(bundle)
    integer(index, 'task index', 0, len(manifest['tasks']) - 1)
    task = manifest['tasks'][index]
    validate_allocation(env, manifest['resources'])
    engine = checked_env_path(env, 'ORCA_EXE')
    scratch = checked_env_path(env, 'HPC_SCRATCH', directory=True)
    result = bundle / 'results' / task['id']
    result.mkdir(parents=True, exist_ok=False)
    native = Path(tempfile.mkdtemp(prefix='pincer-' + task['id'] + '-', dir=scratch))
    for name in ('job.inp', 'geometry.xyz', 'metadata.json'):
        shutil.copyfile(bundle / 'tasks' / task['id'] / name, native / name)
    receipt = dict(schema='pincer_hpc_run_v1', status='running', executed=False,
                   task_id=task['id'], manifest_sha256=sha256(bundle / 'manifest.json'),
                   protocol_sha256=manifest['protocol_sha256'], input_sha256=task['input_sha256'],
                   geometry_sha256=task['geometry_sha256'], engine_binary_sha256=sha256(engine),
                   engine_fingerprint_scope='ORCA driver binary; native output version; not all MPI/shared libraries',
                   engine_path=str(engine), start_utc=datetime.now(timezone.utc).isoformat(),
                   scratch_path=str(native), slurm={k: env[k] for k in
                   ('SLURM_JOB_ID', 'SLURM_NTASKS', 'SLURM_JOB_NUM_NODES', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_NODE')},
                   loaded_modules=env.get('LOADEDMODULES'), fixture=False,
                   physical_kinetic_certificate=False)
    dump(result / 'run.json', receipt)
    env.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    env['PATH'] = str(engine.parent) + os.pathsep + env.get('PATH', '')
    try:
        with (native / 'job.out').open('wb') as out, (native / 'job.err').open('wb') as err:
            proc = subprocess.run([str(engine), 'job.inp'], cwd=native, env=env,
                                  stdout=out, stderr=err, check=False)
        receipt.update(status='completed', executed=True, returncode=proc.returncode)
    except Exception as exc:
        receipt.update(status='runner_failed', failure=type(exc).__name__ + ': ' + str(exc))
    finally:
        # Keep scratch; no destructive cleanup. A hard scheduler kill can leave a
        # running receipt, which is intentionally uncollectable until recovered.
        for p in native.rglob('*'):
            if p.is_symlink():
                raise ValueError('refuse symlink in native calculation directory')
        shutil.copytree(native, result / 'native')
        receipt['native_files'] = {p.relative_to(result / 'native').as_posix(): sha256(p)
                                   for p in sorted((result / 'native').rglob('*')) if p.is_file()}
        receipt['end_utc'] = datetime.now(timezone.utc).isoformat()
        dump(result / 'run.json', receipt)
    return receipt


def verify_result(bundle, task, allow_fixture=False):
    bundle = Path(bundle)
    folder = bundle / 'results' / task['id']
    receipt = load_json(relative_file(bundle, f'results/{task["id"]}/run.json'))
    if receipt.get('schema') != 'pincer_hpc_run_v1' or receipt.get('status') != 'completed' or receipt.get('executed') is not True:
        raise ValueError('missing completed execution status; results are unknown')
    if receipt.get('fixture') is not False and not allow_fixture:
        raise ValueError('test fixtures cannot be collected as physical results')
    if type(receipt.get('returncode')) is not int:
        raise ValueError('missing native exit code')
    if receipt.get('task_id') != task['id'] or receipt.get('manifest_sha256') != sha256(bundle / 'manifest.json'):
        raise ValueError('run / bundle identity mismatch')
    if receipt.get('protocol_sha256') != digest_object(PROTOCOL):
        raise ValueError('run protocol mismatch')
    for key in ('input_sha256', 'geometry_sha256'):
        if receipt.get(key) != task[key]:
            raise ValueError('run source identity mismatch')
    if not isinstance(receipt.get('engine_binary_sha256'), str) or not SHA.fullmatch(receipt['engine_binary_sha256']):
        raise ValueError('missing engine binary fingerprint')
    native = folder / 'native'
    declared = receipt.get('native_files')
    if not isinstance(declared, dict) or not {'job.inp', 'job.out', 'job.err', 'geometry.xyz', 'metadata.json'} <= set(declared):
        raise ValueError('missing native evidence index')
    actual = {p.relative_to(native).as_posix() for p in native.rglob('*') if p.is_file()}
    if set(declared) != actual:
        raise ValueError('native evidence membership mismatch')
    for name, expected in declared.items():
        if sha256(relative_file(bundle, f'results/{task["id"]}/native/{name}')) != expected:
            raise ValueError('native evidence hash mismatch')
    if declared['job.inp'] != task['input_sha256'] or declared['geometry.xyz'] != task['geometry_sha256'] or declared['metadata.json'] != task['metadata_sha256']:
        raise ValueError('native input changed from planned input')
    parsed = parse_orca_output((native / 'job.out').read_text(encoding='utf-8', errors='replace'), task['multiplicity'])
    if parsed['fixture_marker_detected'] and not allow_fixture:
        raise ValueError('software fixture output cannot become physical evidence')
    if receipt['returncode'] != 0:
        parsed['task_converged'] = False
        parsed['eligible_for_vertical_spin_diagnostic'] = False
    return dict(task_id=task['id'], run_sha256=sha256(folder / 'run.json'),
                fixture=receipt.get('fixture'), returncode=receipt['returncode'], **parsed)


def submit(bundle, mode, pilot_review=None, env=None, execute=False):
    env = dict(os.environ if env is None else env)
    bundle = Path(bundle).resolve()
    manifest = verify_bundle(bundle)
    for key in ('HPC_ACCOUNT', 'HPC_PARTITION'):
        if not SAFE_ID.fullmatch(env.get(key, '')):
            raise ValueError(f'{key} requires explicit valid site value')
    for key in ('ORCA_EXE', 'HPC_PYTHON'):
        checked_env_path(env, key)
    checked_env_path(env, 'HPC_SCRATCH', directory=True)
    if mode == 'pilot':
        array = '0-0%1'
    elif mode == 'production':
        if not pilot_review:
            raise ValueError('production requires a human pilot review receipt')
        review = load_json(pilot_review)
        pilot = verify_result(bundle, manifest['tasks'][0])
        if (review.get('approved') is not True or not review.get('reviewer')
                or not review.get('reviewed_utc') or not pilot['task_converged']
                or review.get('manifest_sha256') != sha256(bundle / 'manifest.json')
                or review.get('pilot_run_sha256') != pilot['run_sha256']):
            raise ValueError('pilot execution or human review not bound to this bundle')
        if len(manifest['tasks']) < 2:
            raise ValueError('no tasks remain after pilot')
        array = f'1-{len(manifest["tasks"]) - 1}%{manifest["resources"]["concurrency"]}'
    else:
        raise ValueError('mode must be pilot or production')
    command = ['sbatch', '--parsable', '--account=' + env['HPC_ACCOUNT'],
               '--partition=' + env['HPC_PARTITION'], '--array=' + array, 'array.sbatch']
    if execute:
        process = subprocess.run(command, cwd=bundle, env=env, capture_output=True, text=True, check=True)
        receipt = dict(command=command, stdout=process.stdout, stderr=process.stderr,
                       submitted_utc=datetime.now(timezone.utc).isoformat(),
                       manifest_sha256=sha256(bundle / 'manifest.json'))
        dump(bundle / ('submission-' + mode + '-' + uuid.uuid4().hex + '.json'), receipt)
        return receipt
    return dict(command=command, submitted=False)


def collect(bundle, output):
    bundle, output = Path(bundle).resolve(), Path(output).resolve()
    manifest = verify_bundle(bundle)
    if output.exists():
        raise FileExistsError('immutable collection destination already exists')
    rows = []
    present = []
    for task in manifest['tasks']:
        folder = bundle / 'results' / task['id']
        if not folder.exists():
            rows.append(dict(task_id=task['id'], status='not_returned', physical_kinetic_certificate=False))
            continue
        row = verify_result(bundle, task)
        rows.append(dict(status='collected_diagnostic', **row))
        present.append(folder)
    if not present:
        raise ValueError('no executed result returned; prepared jobs are not calculations')
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.parent / (output.name + '.collecting-' + uuid.uuid4().hex)
    staging.mkdir()
    files = [bundle / 'manifest.json'] + [relative_file(bundle, name) for name in manifest['files']]
    files += [p for folder in present for p in sorted(folder.rglob('*')) if p.is_file()]
    # Return every raw native file, not just accepted rows; never extract a user archive.
    files += sorted(bundle.glob('slurm-*.out')) + sorted(bundle.glob('slurm-*.err'))
    files += sorted(bundle.glob('submission-*.json')) + sorted(bundle.glob('sacct*.txt'))
    for name in ('pilot_review.json', 'site_information.json'):
        if (bundle / name).is_file():
            files.append(bundle / name)
    archive = staging / 'hpc_return.zip'
    archive_index = {}
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zipped:
        for file in files:
            name = file.relative_to(bundle).as_posix()
            file = relative_file(bundle, name)
            before = sha256(file)
            zipped.write(file, name)
            if sha256(file) != before:
                raise ValueError('source changed while archiving')
            archive_index[name] = before
    with zipfile.ZipFile(archive) as zipped:
        if zipped.testzip() is not None or set(zipped.namelist()) != set(archive_index):
            raise ValueError('return archive validation failed')
        for name, expected in archive_index.items():
            h = hashlib.sha256()
            with zipped.open(name) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(chunk)
            if h.hexdigest() != expected:
                raise ValueError('return archive member hash mismatch')
    result = dict(schema='pincer_hpc_collection_v1', rows=rows,
                  manifest_sha256=sha256(bundle / 'manifest.json'),
                  protocol_sha256=manifest['protocol_sha256'], archive_sha256=sha256(archive),
                  archive_files=archive_index, physical_kinetic_certificate=False,
                  limitation='Vertical electronic diagnostics only; no validated minimum, TS, MECP, or rate.')
    dump(staging / 'collection.json', result)
    staging.rename(output)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('verify', 'run-task', 'submit', 'collect'):
        child = commands.add_parser(name)
        child.add_argument('--bundle', type=Path, required=True)
        if name == 'run-task':
            child.add_argument('--index', type=int, required=True)
        if name == 'submit':
            child.add_argument('--mode', choices=('pilot', 'production'), required=True)
            child.add_argument('--pilot-review', type=Path)
            child.add_argument('--execute', action='store_true', help='Actually call sbatch on this site')
        if name == 'collect':
            child.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == 'verify':
        result = verify_bundle(args.bundle)
        result = dict(status='offline_bundle_verified_not_executed', tasks=len(result['tasks']),
                      manifest_sha256=sha256(args.bundle / 'manifest.json'))
    elif args.command == 'run-task':
        result = run_task(args.bundle, args.index)
    elif args.command == 'submit':
        result = submit(args.bundle, args.mode, args.pilot_review, execute=args.execute)
    else:
        result = collect(args.bundle, args.output)
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    if args.command == 'run-task':
        return 0 if result.get('status') == 'completed' and result.get('returncode') == 0 else 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
