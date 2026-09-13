"""Analytic software tests; these surfaces are never catalyst training labels."""
import json
import subprocess

import numpy as np
import pytest

from pincer_catmech.quantum.spin_mecp import (
    MECPConfig, SurfacePair, damped_bfgs, harvey_effective_gradient,
    kkt_step, optimize_mecp, projected_gradient, saturated_penalty, seam_curvature,
)
from pincer_catmech.quantum.spin_provider import Psi4SubprocessProvider, validate_spin_system


def parabola_pair(x):
    z = x.ravel()
    center = np.zeros_like(z)
    center[0] = 1
    a, b = z + center, z - center
    return SurfacePair(float(a @ a / 2), float(b @ b / 2), a.reshape(x.shape),
                       b.reshape(x.shape), "analytic_state_1", "analytic_state_3", True,
                       {"evidence": "synthetic analytic verification only"})


def test_exact_regular_crossing():
    result = optimize_mecp(parabola_pair, np.array([[.8, .6, -.3]]),
                           MECPConfig(gap_tolerance_hartree=1e-9, projected_gradient_tolerance=1e-9))
    assert result.converged, result.to_dict()
    assert abs(result.gap_hartree) < 1e-9
    assert result.projected_gradient_max < 1e-9
    assert np.linalg.norm(result.positions_bohr) < 1e-8
    assert not result.minimum_verified
    json.dumps(result.to_dict(), allow_nan=False)


def test_kkt_step_solves_block_equations():
    b = np.array([[2., .2, 0], [.2, 3., .1], [0., .1, 1.]])
    g, d, gap = np.array([.1, -.4, .3]), np.array([1., 2., -.2]), .7
    step, multiplier = kkt_step(b, g, d, gap)
    assert np.max(np.abs(b @ step + d * multiplier + g)) < 1e-13
    assert abs(d @ step + gap) < 1e-13


def test_projector_and_lagrange_conditions_agree():
    pair = parabola_pair(np.array([[.2, .6, -.3]]))
    tangent, multiplier = projected_gradient(pair)
    assert abs(pair.difference_gradient @ tangent) < 1e-14
    np.testing.assert_allclose(tangent, pair.average_gradient + multiplier * pair.difference_gradient)
    assert np.isclose(np.linalg.norm(harvey_effective_gradient(pair))**2,
                      np.linalg.norm(tangent)**2 + pair.gap**2 * np.linalg.norm(pair.difference_gradient)**2)


def test_penalty_gradient_matches_central_difference():
    x = np.array([[.23, -.12, .4]])
    _, analytic = saturated_penalty(parabola_pair(x), alpha=4, epsilon=.07)
    numeric = np.zeros(x.size)
    for j in range(x.size):
        delta = np.eye(x.size)[j].reshape(x.shape) * 1e-6
        fp = saturated_penalty(parabola_pair(x+delta), 4, .07)[0]
        fm = saturated_penalty(parabola_pair(x-delta), 4, .07)[0]
        numeric[j] = (fp-fm)/2e-6
    np.testing.assert_allclose(analytic, numeric, atol=2e-9, rtol=0)


def test_negative_gap_is_not_accepted_on_inequality_alone():
    result = optimize_mecp(parabola_pair, np.array([[-2., 0, 0]]), MECPConfig(max_iterations=1))
    assert result.gap_hartree < -1e-4
    assert not result.converged


@pytest.mark.parametrize("spin_dependent,state2", [(False, "b"), (True, "a")])
def test_unphysical_or_same_state_provider_rejected(spin_dependent, state2):
    def provider(x):
        return SurfacePair(0, 0, x, x, "a", state2, spin_dependent)
    result = optimize_mecp(provider, np.ones((2, 3)))
    assert not result.converged
    assert result.status == "invalid_surfaces_or_backend_failure"


def test_equal_energies_and_gradients_not_a_crossing_certificate():
    def provider(x):
        return SurfacePair(0, 0, np.zeros_like(x), np.zeros_like(x), "s", "t", True)
    result = optimize_mecp(provider, np.ones((2, 3)))
    assert result.status == "degenerate_surfaces"
    assert not result.converged


def test_parallel_separated_surfaces_singular():
    def provider(x):
        return SurfacePair(-1, 1, x, x, "s", "t", True)
    result = optimize_mecp(provider, np.ones((2, 3)))
    assert result.status == "singular_constraint"
    assert not result.converged


def test_actual_gradient_is_required_to_have_same_shape():
    def provider(x):
        return SurfacePair(0, 1, np.zeros(6), np.zeros(6), "s", "t", True)
    result = optimize_mecp(provider, np.ones((2, 3)))
    assert not result.converged and "match" in result.error


def test_backend_failure_remains_failure():
    def provider(x):
        raise RuntimeError("SCF did not converge")
    result = optimize_mecp(provider, np.ones((2, 3)))
    assert not result.converged
    assert result.energy1_hartree is None
    assert "SCF" in result.error


def test_budget_is_strict_on_number_of_pairs():
    result = optimize_mecp(parabola_pair, np.ones((2, 3)), MECPConfig(max_evaluations=1))
    assert result.status == "budget_exhausted"
    assert result.evaluations == 1


def test_damped_bfgs_handles_negative_secant_curvature():
    b = np.eye(3)
    step, y = np.array([1., 2., 3.]), np.array([-4., -.3, -2.])
    updated = damped_bfgs(b, step, y)
    assert np.linalg.eigvalsh(updated).min() > 0
    np.testing.assert_allclose(updated, updated.T, atol=1e-15)


def radial_pair(x):
    vector = x[1] - x[0]
    r = np.linalg.norm(vector)
    unit = vector / r
    g1, g2 = np.stack((-unit, unit))*(r-1), np.stack((-unit, unit))*(r-3)
    return SurfacePair(.5*(r-1)**2, .5*(r-3)**2, g1, g2, "s", "t", True)


def test_crossing_and_trajectory_covariant_under_translation_and_rotation():
    initial = np.array([[0., 0., 0.], [1.4, .3, -.2]])
    q, _ = np.linalg.qr(np.array([[.4, .7, -.1], [.2, -.4, .8], [.9, .1, .3]]))
    shift = np.array([19., -8., 3.])
    cfg = MECPConfig(gap_tolerance_hartree=1e-10, projected_gradient_tolerance=1e-10)
    a = optimize_mecp(radial_pair, initial, cfg)
    b = optimize_mecp(radial_pair, initial @ q.T + shift, cfg)
    assert a.converged and b.converged
    np.testing.assert_allclose(b.positions_bohr, a.positions_bohr @ q.T + shift, atol=1e-9, rtol=0)
    assert np.isclose(np.linalg.norm(a.positions_bohr[1]-a.positions_bohr[0]), 2.)


def test_curved_seam_converges_to_kkt_point():
    def provider(x):
        z = x.ravel()
        g = z - np.array([.25, .8, 0])
        average = float(g @ g / 2)
        gap = z[0] + .3*z[1]**2 - .1
        d = np.array([1., .6*z[1], 0])
        return SurfacePair(average+gap/2, average-gap/2, (g+d/2).reshape(x.shape),
                           (g-d/2).reshape(x.shape), "s", "t", True)
    result = optimize_mecp(provider, np.array([[.4, .1, .2]]),
                           MECPConfig(max_iterations=100, gap_tolerance_hartree=1e-8,
                                      projected_gradient_tolerance=1e-7))
    assert result.converged, result.to_dict()
    assert abs(result.gap_hartree) < 1e-8
    assert result.projected_gradient_max < 1e-7


def test_curvature_certificate_separates_minimum_from_saddle():
    minimum = seam_curvature(parabola_pair, np.zeros((1, 3)), remove_rigid_motion=False)
    assert minimum["minimum_verified"]
    np.testing.assert_allclose(minimum["eigenvalues_hartree_bohr2"], [1, 1], atol=1e-12)
    away = seam_curvature(parabola_pair, np.ones((1, 3)), remove_rigid_motion=False)
    assert away["positive_curvature"] and not away["minimum_verified"]
    def saddle(x):
        pair = parabola_pair(x)
        pair.energy1_hartree -= x[0, 1]**2
        pair.energy2_hartree -= x[0, 1]**2
        pair.gradient1_hartree_bohr[0, 1] -= 2*x[0, 1]
        pair.gradient2_hartree_bohr[0, 1] -= 2*x[0, 1]
        return pair
    crossing = optimize_mecp(saddle, np.zeros((1, 3)))
    assert crossing.converged and not crossing.minimum_verified
    certificate = seam_curvature(saddle, crossing.positions_bohr, remove_rigid_motion=False)
    assert not certificate["minimum_verified"]
    assert min(certificate["eigenvalues_hartree_bohr2"]) < 0


@pytest.mark.parametrize("charge,multiplicities", [(0, (2, 4)), (0, (1, 1)), (0, (1, 19)), (False, (1, 3))])
def test_spin_electron_parity_and_identities(charge, multiplicities):
    with pytest.raises(ValueError):
        validate_spin_system(["O", "O"], [[0, 0, 0], [0, 0, 2.3]], charge, multiplicities)


def test_even_and_odd_spin_rules():
    validate_spin_system(["Fe", "H", "H"], [[0,0,0], [0,0,3], [0,3,0]], 0, (1,3))
    validate_spin_system(["Co", "H", "H"], [[0,0,0], [0,0,3], [0,3,0]], 0, (2,4))


def test_subprocess_timeout_preserves_launch_evidence(tmp_path, monkeypatch):
    executable = tmp_path / "fake_python"
    executable.write_text("not executed")
    provider = Psi4SubprocessProvider(executable, ["O", "O"], workdir=tmp_path / "native", timeout_seconds=.01)
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], .01)
    monkeypatch.setattr(subprocess, "run", timeout)
    with pytest.raises(TimeoutError):
        provider(np.array([[0,0,0], [0,0,2.3]]))
    records = list((tmp_path / "native").glob("job_*/launch.json"))
    assert len(records) == 1
    assert json.loads(records[0].read_text())["status"] == "timeout"


@pytest.mark.parametrize("settings", [dict(max_iterations=0), dict(wall_seconds=np.nan),
                                      dict(initial_penalty=10, maximum_penalty=1), dict(armijo=1)])
def test_invalid_config(settings):
    with pytest.raises(ValueError):
        MECPConfig(**settings)


def test_spin_cache_requires_raw_hash_and_numeric_agreement(tmp_path):
    import hashlib
    import importlib.util
    from pathlib import Path
    spec=importlib.util.spec_from_file_location("spin_cli",Path(__file__).resolve().parents[1]/"scripts/run_phase3_spin.py")
    cli=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    raw=tmp_path/"response.json"
    state=dict(energy_hartree=-1.2,gradient_hartree_bohr=[[1,2,3]],spin_squared=2.01)
    raw.write_text(json.dumps(dict(states=[state],geometry_sha256="geometry")))
    record=dict(status="converged",signature_sha256="signature",geometry_sha256="geometry",**state,
                native_evidence=[dict(path=str(raw),sha256=hashlib.sha256(raw.read_bytes()).hexdigest())])
    index=tmp_path/"result.json"
    index.write_text(json.dumps(record))
    assert cli.accepted_cache(index,"signature") is not None
    assert cli.accepted_cache(index,"wrong signature") is None
    record["energy_hartree"]=-99
    index.write_text(json.dumps(record))
    assert cli.accepted_cache(index,"signature") is None
    record["energy_hartree"]=-1.2
    index.write_text(json.dumps(record))
    raw.write_text("tampered native output")
    assert cli.accepted_cache(index,"signature") is None


def test_vertical_gaps_require_matched_protocol_and_flag_spin_purity(tmp_path):
    import importlib.util
    from pathlib import Path
    from types import SimpleNamespace
    spec=importlib.util.spec_from_file_location("spin_cli_gaps",Path(__file__).resolve().parents[1]/"scripts/run_phase3_spin.py")
    cli=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    args=SimpleNamespace(method="pbe",basis="sto-3g",chemical_state="active",grid_radial=35,grid_spherical=110,
                         scf_algorithm="soscf",soscf_start_convergence=1e-4,planned_states=2,
                         e_convergence=1e-8,d_convergence=1e-6,
                         output_dir=tmp_path,aggregate_path=tmp_path/"summary.json")
    a=dict(catalyst_id="Fe_example",multiplicity=1,status="converged",energy_hartree=-10.,charge=0,
           geometry_sha256="same_geometry",paired_protocol_sha256="coarse_grid",signature_sha256="a",
           spin_contamination_flag=False,intended_catalyst_identity_pass=True)
    b={**a,"multiplicity":3,"energy_hartree":-9.99,"signature_sha256":"b","paired_protocol_sha256":"fine_grid"}
    summary=cli.summarize(args,[a,b],[])
    assert summary["spin_gaps"][0]["status"]=="missing_or_failed_state"
    b["paired_protocol_sha256"]="coarse_grid"
    summary=cli.summarize(args,[a,b],[])
    gap=summary["spin_gaps"][0]
    assert gap["diagnostic_label_eligible"] and not gap["validated_catalyst_label"]
    assert np.isclose(gap["gap_high_minus_low_hartree"],.01)
    b["spin_contamination_flag"]=True
    summary=cli.summarize(args,[a,b],[])
    assert not summary["spin_gaps"][0]["diagnostic_label_eligible"]
