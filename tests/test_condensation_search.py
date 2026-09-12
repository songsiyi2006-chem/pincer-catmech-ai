"""Geometry/numerical unit cases; these are not computed research observations."""
import numpy as np
import pytest
from ase import Atoms

from pincer_catmech.kinetics.condensation_search import (
    condensation_mode, coordinate_gradient, identity, observed_bonds,
)


def test_covalent_identity_requires_graph_and_carbon_nitrogen_metric():
    atoms = Atoms("HCN", positions=[[-1.06, 0, 0], [0, 0, 0], [1.3, 0, 0]])
    expected = {(0, 1), (1, 2)}
    assert observed_bonds(atoms) == expected
    assert identity(atoms, expected, cn=(1, 2), cn_interval=(1.15, 1.42))["retained"]
    assert not identity(atoms, expected, cn=(1, 2), cn_interval=(1.40, 1.70))["retained"]
    assert not identity(atoms, {(1, 2)}, cn=(1, 2), cn_interval=(1.15, 1.42))["retained"]


def test_reaction_coordinate_gradient_matches_independent_central_differences():
    atoms = Atoms("NHO", positions=[[0, 0, 0], [1.1, 0.2, 0.1], [2.8, 0.4, 0.2]])
    terms = [(0, 1, 1), (2, 1, -1)]
    calculated = coordinate_gradient(atoms, terms)
    finite = np.zeros((3, 3))
    step = 1e-6
    for i in range(3):
        for j in range(3):
            shifted = atoms.copy()
            shifted.positions[i, j] += step
            plus = shifted.get_distance(0, 1)-shifted.get_distance(2, 1)
            shifted.positions[i, j] -= 2*step
            minus = shifted.get_distance(0, 1)-shifted.get_distance(2, 1)
            finite[i, j] = (plus-minus)/(2*step)
    np.testing.assert_allclose(calculated, finite, atol=5e-10)
    np.testing.assert_allclose(calculated.sum(axis=0), 0, atol=1e-14)


def test_mode_sign_is_arbitrary_but_transverse_motion_is_rejected():
    atoms = Atoms("NHO", positions=[[0, 0, 0], [1.1, 0, 0], [2.8, 0, 0]])
    spec = {"transfer": {"terms": [(0, 1, 1), (2, 1, -1)], "minimum_overlap": 0.15}}
    longitudinal = np.array([[-1., 0, 0], [2., 0, 0], [-1., 0, 0]])
    assert condensation_mode(atoms, longitudinal, spec)["accepted"]
    assert condensation_mode(atoms, -longitudinal, spec)["accepted"]
    transverse = longitudinal[:, [1, 0, 2]]
    assert not condensation_mode(atoms, transverse, spec)["accepted"]


@pytest.mark.parametrize("mode", [np.zeros((3, 3)), np.full((3, 3), np.nan)])
def test_invalid_mode_cannot_be_certified(mode):
    atoms = Atoms("NHO", positions=[[0, 0, 0], [1.1, 0, 0], [2.8, 0, 0]])
    with pytest.raises(ValueError):
        condensation_mode(atoms, mode, {"transfer": {"terms": [(0, 1, 1)], "minimum_overlap": 0.1}})
