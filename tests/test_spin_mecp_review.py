"""Adversarial mathematical/provenance tests, never electronic-structure data."""
from copy import deepcopy

import numpy as np
import pytest

from pincer_catmech.quantum.spin_mecp import (
    SurfacePair, harvey_effective_gradient, optimize_mecp, projected_gradient,
    seam_curvature,
)


def quadratic_pair(x):
    center = np.zeros_like(x)
    center.flat[0] = 1
    a, b = x + center, x - center
    return SurfacePair(float(np.sum(a*a)/2), float(np.sum(b*b)/2), a, b,
                       "multiplicity_1", "multiplicity_3", True)


def with_metadata(x, metadata):
    pair = quadratic_pair(x)
    pair.metadata = deepcopy(metadata)
    return pair


def spin_metadata(s2=2.01):
    return {"states": [
        {"multiplicity": 1, "spin_squared": 0., "expected_spin_squared": 0., "scf_converged": True},
        {"multiplicity": 3, "spin_squared": s2, "expected_spin_squared": 2., "scf_converged": True},
    ]}


@pytest.mark.parametrize("s2", [2.6614, 1.8, np.nan, np.inf, -1., [2.]])
def test_contaminated_or_malformed_supplied_spin_is_not_accepted(s2):
    result = optimize_mecp(lambda x: with_metadata(x, spin_metadata(s2)), np.zeros((1, 3)))
    assert result.status == "invalid_surfaces_or_backend_failure"
    assert not result.converged
    with pytest.raises(ValueError):
        seam_curvature(lambda x: with_metadata(x, spin_metadata(s2)), np.zeros((1, 3)),
                       remove_rigid_motion=False)


def test_spin_screen_is_eligibility_only_and_metadata_free_toys_still_work():
    for provider in (quadratic_pair, lambda x: with_metadata(x, spin_metadata(2.099))):
        result = optimize_mecp(provider, np.zeros((1, 3)))
        assert result.converged and not result.minimum_verified
    metadata = spin_metadata()
    del metadata["states"][1]["expected_spin_squared"]
    assert optimize_mecp(lambda x: with_metadata(x, metadata), np.zeros((1, 3))).converged


@pytest.mark.parametrize("quality", [
    {"scf_converged": False}, {"scf_converged": "false"},
    {"spin_contamination_flag": True}, {"intended_catalyst_identity_pass": False},
    {"diagnostic_label_eligible": False}, {"status": "failed"}, {"quality_pass": False},
    {"gradient_kind": "finite_difference"},
])
@pytest.mark.parametrize("location", ["pair", "state"])
def test_supplied_failed_quality_never_certifies_a_crossing(quality, location):
    metadata = quality if location == "pair" else {"states": [{}, quality]}
    result = optimize_mecp(lambda x: with_metadata(x, metadata), np.zeros((1, 3)))
    assert not result.converged and result.error


@pytest.mark.parametrize("change", [
    {"multiplicity": 5}, {"expected_spin_squared": 2.01},
    {"energy_hartree": 99.}, {"gradient_hartree_bohr": [[99., 0., 0.]]},
])
def test_state_metadata_must_match_returned_state(change):
    metadata = spin_metadata()
    metadata["states"][1].update(change)
    result = optimize_mecp(lambda x: with_metadata(x, metadata), np.zeros((1, 3)))
    assert not result.converged and result.status == "invalid_surfaces_or_backend_failure"


@pytest.mark.parametrize("metadata", [
    {"coordinate_units": "angstrom"}, {"energy_units": "eV"},
    {"gradient_units": "eV/angstrom"}, {"positions_bohr": [[1., 0., 0.]]},
    {"states": [{"charge": 0}, {"charge": 1}]},
    {"method": "pbe", "states": [{"method": "hf"}, {}]},
    {"states": [{}]}, {"states": [None, {}]},
])
def test_mismatched_geometry_units_and_pair_provenance_are_rejected(metadata):
    result = optimize_mecp(lambda x: with_metadata(x, metadata), np.zeros((1, 3)))
    assert result.status == "invalid_surfaces_or_backend_failure"


@pytest.mark.parametrize("mutation", ["basis", "state", "metadata_disappears", "reference"])
def test_protocol_must_stay_fixed_during_search_and_curvature(mutation):
    def make_provider():
        calls = 0
        def provider(x):
            nonlocal calls
            calls += 1
            metadata = {"basis": "test", "states": [{"reference": "rks"}, {"reference": "uks"}]}
            if calls > 1:
                if mutation == "basis":
                    metadata["basis"] = "other"
                elif mutation == "metadata_disappears":
                    del metadata["states"]
                elif mutation == "reference":
                    metadata["states"][1]["reference"] = "uhf"
            pair = with_metadata(x, metadata)
            if calls > 1 and mutation == "state":
                pair.state2 = "another_state"
            return pair
        return provider
    result = optimize_mecp(make_provider(), np.ones((1, 3)))
    assert result.status == "invalid_surfaces_or_backend_failure"
    with pytest.raises(ValueError):
        seam_curvature(make_provider(), np.zeros((1, 3)), remove_rigid_motion=False)


def test_symmetrization_cannot_certify_a_nonconservative_gradient_field():
    def provider(x):
        pair = quadratic_pair(x)
        skew = np.array([[0., 20.*x[0, 2], -20.*x[0, 1]]])
        pair.gradient1_hartree_bohr += skew
        pair.gradient2_hartree_bohr += skew
        return pair
    result = seam_curvature(provider, np.zeros((1, 3)), remove_rigid_motion=False)
    assert result["positive_curvature"]  # The symmetric part deceptively looks ideal.
    np.testing.assert_allclose(result["eigenvalues_hartree_bohr2"], [1., 1.])
    assert result["antisymmetry_max"] == pytest.approx(40.)
    assert not result["hessian_symmetry_pass"] and not result["minimum_verified"]


def test_curved_constraint_requires_the_lagrangian_not_just_energy_hessian():
    def provider(x):
        a, b, c = x.ravel()
        average = -2*a - .5*b*b + .5*c*c
        gap = a + .5*b*b
        g, d = np.array([-2., -b, c]), np.array([1., b, 0.])
        return SurfacePair(average+gap/2, average-gap/2, (g+d/2).reshape(x.shape),
                           (g-d/2).reshape(x.shape), "a", "b", True)
    pair = provider(np.zeros((1, 3)))
    tangent, multiplier = projected_gradient(pair)
    assert multiplier == 2 and np.max(abs(tangent)) == 0
    result = seam_curvature(provider, np.zeros((1, 3)), remove_rigid_motion=False)
    # Along a=-b^2/2, A=b^2/2+c^2/2, although A_yy=-1 in Cartesian space.
    np.testing.assert_allclose(result["eigenvalues_hartree_bohr2"], [1., 1.], atol=1e-12)
    assert result["minimum_verified"]


def test_constraint_is_not_lost_against_large_rotation_vectors():
    origin = np.array([[0., 0., 0.], [1e8, 0., 0.], [0., 1e8, 0.]])
    d = np.array([[-1e-9, 0., 0.], [1e-9, 0., 0.], [0., 0., 0.]])
    def provider(x):
        delta = x-origin
        a, c = float(np.sum(delta*delta)/2), float(np.sum(d*delta))
        return SurfacePair(a+c/2, a-c/2, delta+d/2, delta-d/2, "a", "b", True)
    result = seam_curvature(provider, origin)
    assert result["tangent_dimension"] == 2  # 9 coordinates - 6 rigid modes - 1 constraint.


@pytest.mark.parametrize("settings", [
    {"antisymmetry_tolerance": -1.}, {"antisymmetry_tolerance": np.inf},
    {"curvature_tolerance": -1.}, {"curvature_tolerance": np.nan},
    {"gap_tolerance_hartree": np.inf}, {"projected_gradient_tolerance": -1.},
    {"displacement_bohr": np.nan},
])
def test_invalid_certificate_tolerances_fail_before_backend_call(settings):
    def forbidden(x):
        pytest.fail("Invalid settings must be rejected before calling a backend")
    with pytest.raises(ValueError):
        seam_curvature(forbidden, np.zeros((1, 3)), **settings)


@pytest.mark.parametrize("x", [[], [0, 0, 0], [[np.nan, 0, 0]], [[0, 0]]])
def test_invalid_curvature_geometry_is_rejected(x):
    with pytest.raises(ValueError):
        seam_curvature(quadratic_pair, x)


@pytest.mark.parametrize("value", [np.nan, np.inf, 0., -1.])
def test_projection_scales_must_be_positive_finite(value):
    pair = quadratic_pair(np.zeros((1, 3)))
    with pytest.raises(ValueError):
        projected_gradient(pair, value)
    with pytest.raises(ValueError):
        harvey_effective_gradient(pair, value)
