"""Analytic geometry/volume checks; these fixtures are synthetic, not QM data."""

import numpy as np
import pytest
from ase import Atoms

from pincer_catmech.features import BONDI_RADII, bite_angle, buried_volume


@pytest.mark.parametrize("angle", [0, 30, 60, 90, 120, 180])
def test_analytic_bite_angles(angle):
    theta = np.radians(angle)
    atoms = Atoms("RuPP", positions=[[0, 0, 0], [2, 0, 0], [3 * np.cos(theta), 3 * np.sin(theta), 0]])
    assert bite_angle(atoms, 0, 1, 2) == pytest.approx(angle, abs=1e-6)
    assert bite_angle(atoms.positions, 0, 1, 2, symbols=["Ru", "P", "P"]) == pytest.approx(angle, abs=1e-6)


def test_bite_angle_translation_and_rotation():
    xyz = np.array([[0, 0, 0], [1, 0, 1], [0, 1, 0]], dtype=float)
    rotation = np.array([[0, -1, 0], [0, 0, 1], [-1, 0, 0]])
    expected = bite_angle(xyz, 0, 1, 2, symbols=["Mn", "N", "P"])
    assert bite_angle(xyz @ rotation + [80, -30, 4], 0, 1, 2, symbols=["Mn", "N", "P"]) == pytest.approx(expected)


def test_exact_sphere_volume_and_chunk_seed_invariance():
    # A synthetic carbon sphere at the sampling center: exact p=(r/R)^3.
    atoms = Atoms("RuC", positions=np.zeros((2, 3)))
    options = dict(ligand_indices=[0, 1], radius_scale=1, sphere_radius=3.5, n_samples=120_000, seed=72)
    result = buried_volume(atoms, 0, chunk_size=1237, **options)
    exact = 100 * (BONDI_RADII["C"] / 3.5) ** 3
    assert abs(result.percent_buried_volume - exact) < 4 * result.standard_error_percent
    assert result.confidence_interval_percent[0] < exact < result.confidence_interval_percent[1]
    assert result.atom_indices == (1,)
    assert result.effective_radii == (1.7,)
    assert result.n_buried == buried_volume(atoms, 0, chunk_size=3000, **options).n_buried


def test_partial_intersection_two_unit_spheres():
    # Two unit spheres separated by d=1 have lens volume 5*pi/12;
    # dividing by the sampling sphere 4*pi/3 gives occupancy 5/16.
    atoms = Atoms("FeC", positions=[[0, 0, 0], [1, 0, 0]])
    result = buried_volume(atoms, 0, ligand_indices=[1], sphere_radius=1, radius_scale=1,
                           radii_overrides={"C": 1}, n_samples=150_000, seed=53)
    assert abs(result.percent_buried_volume - 31.25) < 4 * result.standard_error_percent


def test_overlapping_atomic_spheres_count_union_once():
    atoms = Atoms("CoCC", positions=[[0, 0, 0], [2, 0, 0], [2, 0, 0]])
    single = buried_volume(atoms, 0, ligand_indices=[1], seed=33, n_samples=20_000)
    overlap = buried_volume(atoms, 0, ligand_indices=[1, 2], seed=33, n_samples=20_000)
    assert single.n_buried == overlap.n_buried
    atoms.translate([11, -6, 2])
    translated = buried_volume(atoms, 0, ligand_indices=[1, 2], seed=33, n_samples=20_000)
    assert translated.n_buried == overlap.n_buried


def test_endpoint_intervals_and_explicit_selection():
    atoms = Atoms("RuCH", positions=[[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    empty = buried_volume(atoms, 0, ligand_indices=[0, 2], n_samples=1000)
    assert empty.atom_indices == ()
    assert empty.percent_buried_volume == empty.standard_error_percent == 0
    assert empty.confidence_interval_percent[0] == 0 < empty.confidence_interval_percent[1]
    full = buried_volume(atoms, 0, ligand_indices=[1], radii_overrides={"C": 10}, n_samples=1000)
    assert full.percent_buried_volume == 100
    assert full.confidence_interval_percent[0] < full.confidence_interval_percent[1] == 100
    hydrogen = buried_volume(atoms, 0, ligand_indices=[2], include_hydrogens=True, n_samples=1000)
    assert hydrogen.percent_buried_volume > 0
    assert hydrogen.effective_radii == pytest.approx((1.2 * 1.17,))


def test_unsupported_ligand_element_requires_override():
    atoms = Atoms("RuMn", positions=[[0, 0, 0], [2, 0, 0]])
    with pytest.raises(ValueError, match="radii_overrides"):
        buried_volume(atoms, 0, ligand_indices=[1])
    result = buried_volume(atoms, 0, ligand_indices=[1], radii_overrides={"Mn": 2.1}, n_samples=100)
    assert result.effective_radii == pytest.approx((2.1 * 1.17,))


def test_xyz_path_and_array_equivalence(tmp_path):
    path = tmp_path / "complex.xyz"
    path.write_text("3\nsynthetic Ru-PNP donor geometry\nRu 0 0 0\nP 2 0 0\nP 0 2 0\n", encoding="utf-8")
    assert bite_angle(path, 0, 1, 2) == pytest.approx(90)
    from_file = buried_volume(path, 0, ligand_indices=[1, 2], n_samples=500)
    from_array = buried_volume([[0, 0, 0], [2, 0, 0], [0, 2, 0]], 0, symbols=["Ru", "P", "P"], ligand_indices=[1, 2], n_samples=500)
    assert from_file == from_array


@pytest.mark.parametrize("text", ["", "no\ncomment\n", "2\ncomment\nC 0 0 0\n", "1\ncomment\nC 0 0\n",
                                  "1\ncomment\nC nan 0 0\n", "1\ncomment\nQq 0 0 0\n",
                                  "1\ncomment\nC 0 0 0\n1\nsecond\nC 0 0 0\n"])
def test_malformed_xyz_rejected(tmp_path, text):
    path = tmp_path / "invalid.xyz"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        buried_volume(path, 0, ligand_indices=[])


@pytest.mark.parametrize("kwargs", [
    {"sphere_radius": 0}, {"sphere_radius": np.inf}, {"radius_scale": -1},
    {"n_samples": 0}, {"n_samples": 2.5}, {"n_samples": True}, {"chunk_size": -1},
    {"seed": -3}, {"confidence_level": 1}, {"confidence_level": np.nan},
    {"ligand_indices": [1, 1]}, {"ligand_indices": [-1]}, {"ligand_indices": [2]},
    {"ligand_indices": [True]}, {"include_hydrogens": "no"}, {"coordinate_units": "bohr"},
    {"radii_overrides": {"C": 0}}, {"radii_overrides": {"Ru": np.inf}},
])
def test_invalid_sampling_arguments(kwargs):
    options = dict(ligand_indices=[1], n_samples=10)
    options.update(kwargs)
    with pytest.raises(ValueError):
        buried_volume(Atoms("RuC", positions=[[0, 0, 0], [2, 0, 0]]), 0, **options)


def test_invalid_geometry_and_indices():
    atoms = Atoms("RuPP", positions=[[0, 0, 0], [0, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError, match="nonzero"):
        bite_angle(atoms, 0, 1, 2)
    with pytest.raises(ValueError, match="distinct"):
        bite_angle(atoms, 0, 2, 2)
    with pytest.raises(ValueError, match="explicit symbols"):
        bite_angle(np.zeros((3, 3)), 0, 1, 2)
    with pytest.raises(ValueError, match="valid element"):
        buried_volume(np.zeros((2, 3)), 0, symbols=["Ru"], ligand_indices=[])
    atoms.pbc = True
    with pytest.raises(ValueError, match="periodic"):
        bite_angle(atoms, 0, 1, 2)
