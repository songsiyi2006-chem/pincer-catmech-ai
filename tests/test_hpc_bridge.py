"""All native-looking outputs here are SOFTWARE FIXTURES, never chemical data."""
from __future__ import annotations
import json
import shutil
from pathlib import Path
import subprocess
import sys
import zipfile

import pytest

from pincer_catmech.hpc import bridge as h


FIXTURE_OUTPUT = '''SOFTWARE FIXTURE -- NOT A QUANTUM CALCULATION
                         Program Version 6.1.1
        **************** SCF CONVERGED AFTER 17 CYCLES ****************
Expectation value of <S**2> : 0.000000001
FINAL SINGLE POINT ENERGY     -1.123456789
                     ****ORCA TERMINATED NORMALLY****
'''


@pytest.fixture
def source(tmp_path):
    src = tmp_path / 'source'
    src.mkdir()
    (src / 'molecule.xyz').write_text('2\nSOFTWARE FIXTURE\nH 0 0 0\nH 0 0 .74\n')
    h.dump(src / 'metadata.json', {'charge': 0, 'formula': 'H2', 'fixture': True})
    selection = src / 'selection.json'
    h.dump(selection, {'schema': 'pincer_hpc_selection_v1', 'records': [
        {'id': 'fixture_h2', 'xyz': 'molecule.xyz', 'metadata': 'metadata.json', 'charge': 0, 'multiplicity': 1}]})
    return src, selection


@pytest.fixture
def resource():
    return h.resources(4, 8192, 1000, '01:00:00', 2)


@pytest.fixture
def bundle(tmp_path, source, resource):
    src, selection = source
    out = tmp_path / 'bundle'
    h.prepare(selection, src, out, resource)
    return out


def fake_result(bundle, fixture=True, output=FIXTURE_OUTPUT, exitcode=0):
    manifest = h.verify_bundle(bundle)
    task = manifest['tasks'][0]
    folder = bundle / 'results' / task['id']
    native = folder / 'native'
    shutil.copytree(bundle / 'tasks' / task['id'], native)
    (native / 'job.out').write_text(output)
    (native / 'job.err').write_text('')
    receipt = dict(schema='pincer_hpc_run_v1', status='completed', executed=True,
                   task_id=task['id'], manifest_sha256=h.sha256(bundle / 'manifest.json'),
                   protocol_sha256=manifest['protocol_sha256'], input_sha256=task['input_sha256'],
                   geometry_sha256=task['geometry_sha256'], engine_binary_sha256='a' * 64,
                   fixture=fixture, returncode=exitcode,
                   native_files={p.name: h.sha256(p) for p in native.iterdir()})
    h.dump(folder / 'run.json', receipt)
    return task, folder


def test_bundle_portable_and_exact_bytes(bundle, source):
    manifest = h.verify_bundle(bundle)
    assert manifest['status'] == 'prepared_not_submitted'
    assert manifest['tasks'][0]['geometry_sha256'] == h.sha256(source[0] / 'molecule.xyz')
    assert not manifest['protocol']['physical_kinetic_certificate']
    assert '%pal nprocs 4 end' in (bundle / 'tasks/fixture_h2/job.inp').read_text()
    assert 'CPCM(toluene)' in (bundle / 'tasks/fixture_h2/job.inp').read_text()
    moved = bundle.with_name('portable_copy')
    shutil.copytree(bundle, moved)
    process = subprocess.run([sys.executable, str(moved / 'bridge.py'), 'verify', '--bundle', str(moved)],
                             capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    assert 'offline_bundle_verified_not_executed' in process.stdout


def test_destination_never_overwritten(bundle, source, resource):
    before = h.sha256(bundle / 'manifest.json')
    with pytest.raises(FileExistsError):
        h.prepare(source[1], source[0], bundle, resource)
    assert h.sha256(bundle / 'manifest.json') == before


@pytest.mark.parametrize('name', ['../secret', '/etc/passwd', 'C:/secret', 'a\\b', './x', 'a/../b'])
def test_path_traversal_rejected(source, name):
    with pytest.raises(ValueError):
        h.relative_file(source[0], name)


@pytest.mark.parametrize('charge,mult', [(0, 2), (True, 1), (0, 0), (0, 31), (50, 1)])
def test_charge_spin_rejection(charge, mult):
    with pytest.raises(ValueError):
        h.validate_electrons([('H', 0, 0, 0), ('H', 0, 0, .74)], charge, mult)


@pytest.mark.parametrize('kwargs', [dict(nprocs=10), dict(maxcore_mb=2000), dict(walltime='0;exit'),
                                   dict(memory_mb=-1), dict(concurrency=0), dict(nprocs=True)])
def test_resource_oversubscribe_or_invalid_rejected(resource, kwargs):
    r = {k: resource[k] for k in ('nprocs', 'memory_mb', 'maxcore_mb', 'walltime', 'concurrency')}
    r.update(kwargs)
    with pytest.raises(ValueError):
        h.resources(**r)


@pytest.mark.parametrize('content', ['1\nfixture\nH nan 0 0\n', '1\nfixture\nH inf 0 0\n',
                                    '2\nfixture\nH 0 0 0\n', '2\nfixture\nH 0 0 0\nH 0 0 .01\n'])
def test_invalid_xyz(source, content):
    p = source[0] / 'bad.xyz'
    p.write_text(content)
    with pytest.raises(ValueError):
        h.read_xyz(p)


def test_metadata_mismatch(source, resource, tmp_path):
    h.dump(source[0] / 'metadata.json', {'charge': 0, 'formula': 'H3'})
    with pytest.raises(ValueError, match='composition'):
        h.prepare(source[1], source[0], tmp_path / 'badbundle', resource)


def test_nonfinite_duplicate_json(source):
    p = source[0] / 'bad.json'
    for value in ['{"a":NaN}', '{"a":1,"a":2}']:
        p.write_text(value)
        with pytest.raises(ValueError):
            h.load_json(p)


def test_bundle_tampering(bundle):
    (bundle / 'tasks/fixture_h2/job.inp').write_text('! fabricated\n')
    with pytest.raises(ValueError, match='hash mismatch'):
        h.verify_bundle(bundle)


def test_parser_conservative_fixture():
    result = h.parse_orca_output(FIXTURE_OUTPUT, 1)
    assert result['task_converged']
    assert result['eligible_for_vertical_spin_diagnostic']
    assert not result['certified_transition_state']
    assert not result['certified_minimum']
    assert result['imaginary_frequency_count'] is None


@pytest.mark.parametrize('removed', ['Program Version 6.1.1', 'SCF CONVERGED AFTER 17 CYCLES',
                                    'FINAL SINGLE POINT ENERGY     -1.123456789', 'ORCA TERMINATED NORMALLY'])
def test_missing_native_status_unknown(removed):
    result = h.parse_orca_output(FIXTURE_OUTPUT.replace(removed, ''), 1)
    assert not result['task_converged']


def test_missing_s2_unknown():
    result = h.parse_orca_output(FIXTURE_OUTPUT.replace('Expectation value of <S**2> : 0.000000001', ''), 1)
    assert result['task_converged']
    assert result['spin_quality'] == 'unknown'
    assert not result['eligible_for_vertical_spin_diagnostic']


def test_contaminated_spin_cannot_be_gap_label():
    result = h.parse_orca_output(FIXTURE_OUTPUT.replace('0.000000001', '2.66144576'), 3)
    assert result['spin_quality'] == 'reject_contaminated'
    assert not result['eligible_for_vertical_spin_diagnostic']


@pytest.mark.parametrize('number', ['nan', 'inf', '-Infinity'])
def test_nonfinite_energy_rejected(number):
    with pytest.raises(ValueError, match='nonfinite'):
        h.parse_orca_output(FIXTURE_OUTPUT.replace('-1.123456789', number), 1)


def test_bad_order_failure_and_multiple_jobs_rejected():
    for text in [FIXTURE_OUTPUT + 'SCF NOT CONVERGED', FIXTURE_OUTPUT + FIXTURE_OUTPUT,
                 FIXTURE_OUTPUT.replace('6.1.1', '5.0.4'),
                 FIXTURE_OUTPUT.replace('FINAL SINGLE POINT ENERGY     -1.123456789', '') +
                 '\nFINAL SINGLE POINT ENERGY -1.123456789\n']:
        assert not h.parse_orca_output(text, 1)['task_converged']


def test_software_fixture_never_collected_as_research(bundle, tmp_path):
    task, _ = fake_result(bundle)
    assert h.verify_result(bundle, task, allow_fixture=True)['fixture']
    with pytest.raises(ValueError, match='fixtures'):
        h.collect(bundle, tmp_path / 'out')


def test_missing_execution_status_rejected(bundle):
    task, folder = fake_result(bundle)
    rec = h.load_json(folder / 'run.json')
    del rec['executed']
    h.dump(folder / 'run.json', rec)
    with pytest.raises(ValueError, match='execution status'):
        h.verify_result(bundle, task, allow_fixture=True)


def test_native_hash_tampering_rejected(bundle):
    task, folder = fake_result(bundle)
    (folder / 'native/job.out').write_text('changed')
    with pytest.raises(ValueError, match='hash mismatch'):
        h.verify_result(bundle, task, allow_fixture=True)


def test_unknown_native_file_rejected(bundle):
    task, folder = fake_result(bundle)
    (folder / 'native/unrecorded.out').write_text('untracked')
    with pytest.raises(ValueError, match='membership'):
        h.verify_result(bundle, task, allow_fixture=True)


def test_native_nonzero_exit_preserved(bundle):
    task, _ = fake_result(bundle, exitcode=1)
    result = h.verify_result(bundle, task, allow_fixture=True)
    assert not result['task_converged']
    assert not result['eligible_for_vertical_spin_diagnostic']


def test_empty_results_not_called_calculated(bundle, tmp_path):
    with pytest.raises(ValueError, match='no executed result'):
        h.collect(bundle, tmp_path / 'out')


def test_site_values_required(bundle):
    with pytest.raises(ValueError, match='HPC_ACCOUNT'):
        h.submit(bundle, 'pilot', env={})


def test_site_command_no_shell_and_pilot_gate(bundle, tmp_path):
    env = dict(HPC_ACCOUNT='lab', HPC_PARTITION='cpu', ORCA_EXE=sys.executable,
               HPC_PYTHON=sys.executable, HPC_SCRATCH=str(tmp_path))
    result = h.submit(bundle, 'pilot', env=env)
    assert '--array=0-0%1' in result['command']
    assert not result['submitted']
    with pytest.raises(ValueError, match='human pilot review'):
        h.submit(bundle, 'production', env=env)
    env['HPC_ACCOUNT'] = 'x; echo secret'
    with pytest.raises(ValueError):
        h.submit(bundle, 'pilot', env=env)


def test_allocation_limits(resource):
    env = dict(SLURM_JOB_ID='1', SLURM_NTASKS='4', SLURM_JOB_NUM_NODES='1',
               SLURM_CPUS_PER_TASK='1', SLURM_MEM_PER_NODE='8192')
    h.validate_allocation(env, resource)
    for key, value in [('SLURM_NTASKS', '2'), ('SLURM_JOB_NUM_NODES', '2'), ('SLURM_MEM_PER_NODE', '1000')]:
        changed = dict(env, **{key: value})
        with pytest.raises(ValueError):
            h.validate_allocation(changed, resource)
    with pytest.raises(ValueError, match='missing scheduler'):
        h.validate_allocation({}, resource)


def test_collect_archive_bytes_without_unlabelled_fixture(bundle, tmp_path, monkeypatch):
    """Mock only the acceptance boundary, retaining fixture=true in every row."""
    task, _ = fake_result(bundle)
    original = h.verify_result
    monkeypatch.setattr(h, 'verify_result', lambda b, t: original(b, t, allow_fixture=True))
    result = h.collect(bundle, tmp_path / 'collected')
    assert result['rows'][0]['fixture'] is True
    assert result['physical_kinetic_certificate'] is False
    archive = tmp_path / 'collected/hpc_return.zip'
    assert h.sha256(archive) == result['archive_sha256']
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist()) == set(result['archive_files'])
        assert b'SOFTWARE FIXTURE' in zipped.read('results/fixture_h2/native/job.out')


def test_crystal_origin_requires_provenance(source, resource, tmp_path):
    selection = h.load_json(source[1])
    selection['records'][0]['structure_origin'] = 'public_crystal_derived_candidate'
    h.dump(source[1], selection)
    with pytest.raises(ValueError, match='source and derivation'):
        h.prepare(source[1], source[0], tmp_path / 'crystal_bundle', resource)


def test_crystal_candidate_never_automatically_certified(source, resource, tmp_path):
    selection = h.load_json(source[1])
    selection['records'][0]['structure_origin'] = 'public_crystal_derived_candidate'
    h.dump(source[1], selection)
    h.dump(source[0] / 'metadata.json', {'charge': 0, 'formula': 'H2', 'fixture': True,
           'provenance': {'source_url': 'https://example.invalid/fixture',
                          'source_identifier': 'SOFTWARE FIXTURE',
                          'structure_derivation': 'Software unit test; no real crystal'}})
    output = tmp_path / 'crystal_bundle'
    manifest = h.prepare(source[1], source[0], output, resource)
    task = manifest['tasks'][0]
    assert task['structure_origin'] == 'public_crystal_derived_candidate'
    assert task['source_provenance']['source_identifier'] == 'SOFTWARE FIXTURE'
    assert not task['ground_spin_verified']
    assert not task['solution_geometry_verified']
    assert 'public_crystal_derived_candidate' in (output / 'README.md').read_text()


def test_origin_tampering_rejected(bundle):
    manifest = h.load_json(bundle / 'manifest.json')
    manifest['tasks'][0]['experimental_identity_verified'] = True
    h.dump(bundle / 'manifest.json', manifest)
    with pytest.raises(ValueError, match='evidence boundary'):
        h.verify_bundle(bundle)


def test_relabelled_fixture_output_still_rejected(bundle):
    task, _ = fake_result(bundle, fixture=False)
    with pytest.raises(ValueError, match='software fixture output'):
        h.verify_result(bundle, task)


def test_runner_preserves_native_outputs_and_refuses_overwrite(bundle, tmp_path, monkeypatch):
    env = dict(ORCA_EXE=sys.executable, HPC_SCRATCH=str(tmp_path), SLURM_JOB_ID='123',
               SLURM_NTASKS='4', SLURM_JOB_NUM_NODES='1', SLURM_CPUS_PER_TASK='1',
               SLURM_MEM_PER_NODE='8192')
    def fake_native(command, cwd, env, stdout, stderr, check):
        assert command == [str(Path(sys.executable).resolve()), 'job.inp']
        assert env['OMP_NUM_THREADS'] == '1'
        stdout.write(FIXTURE_OUTPUT.encode())
        (cwd / 'job.gbw').write_bytes(b'SOFTWARE FIXTURE WAVEFUNCTION')
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(h.subprocess, 'run', fake_native)
    receipt = h.run_task(bundle, 0, env)
    assert receipt['status'] == 'completed'
    assert receipt['executed']
    assert 'job.gbw' in receipt['native_files']
    assert Path(receipt['scratch_path']).exists()
    with pytest.raises(FileExistsError):
        h.run_task(bundle, 0, env)
    task = h.verify_bundle(bundle)['tasks'][0]
    with pytest.raises(ValueError, match='software fixture output'):
        h.verify_result(bundle, task)


def test_runner_launch_failure_does_not_claim_execution(bundle, tmp_path, monkeypatch):
    env = dict(ORCA_EXE=sys.executable, HPC_SCRATCH=str(tmp_path), SLURM_JOB_ID='123',
               SLURM_NTASKS='4', SLURM_JOB_NUM_NODES='1', SLURM_CPUS_PER_TASK='1',
               SLURM_MEM_PER_NODE='8192')
    def fail(*args, **kwargs):
        raise OSError('SOFTWARE FIXTURE launch failure')
    monkeypatch.setattr(h.subprocess, 'run', fail)
    receipt = h.run_task(bundle, 0, env)
    assert receipt['status'] == 'runner_failed'
    assert not receipt['executed']
    assert 'job.inp' in receipt['native_files']


def test_production_skips_pilot_and_binds_review(source, resource, tmp_path, monkeypatch):
    selection = h.load_json(source[1])
    selection['records'].append(dict(selection['records'][0], id='fixture_h2_second'))
    h.dump(source[1], selection)
    bundle = tmp_path / 'two_task_bundle'
    h.prepare(source[1], source[0], bundle, resource)
    task, folder = fake_result(bundle)
    original = h.verify_result
    # Explicit mock boundary only; production calls the unmodified rejecting verifier.
    monkeypatch.setattr(h, 'verify_result', lambda b, t: original(b, t, allow_fixture=True))
    review = tmp_path / 'review.json'
    h.dump(review, dict(approved=True, reviewer='SOFTWARE FIXTURE', reviewed_utc='2026-09-16T00:00:00Z',
                        manifest_sha256=h.sha256(bundle / 'manifest.json'),
                        pilot_run_sha256=h.sha256(folder / 'run.json')))
    env = dict(HPC_ACCOUNT='lab', HPC_PARTITION='cpu', ORCA_EXE=sys.executable,
               HPC_PYTHON=sys.executable, HPC_SCRATCH=str(tmp_path))
    result = h.submit(bundle, 'production', pilot_review=review, env=env)
    assert '--array=1-1%2' in result['command']
    wrong = h.load_json(review)
    wrong['pilot_run_sha256'] = '0' * 64
    h.dump(review, wrong)
    with pytest.raises(ValueError, match='pilot execution or human review'):
        h.submit(bundle, 'production', pilot_review=review, env=env)
