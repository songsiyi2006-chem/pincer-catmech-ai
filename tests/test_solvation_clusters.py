"""Chemical topology, directional docking and balanced energy identities."""
from itertools import combinations
import hashlib
import json
from pathlib import Path
import zipfile

from ase import Atoms
from ase.io import read
import numpy as np
import pytest

from pincer_catmech.generators.solvation_clusters import (
    association_energy, build_solvation_cluster, connectivity_check,
    hydrogen_bond_metrics, many_body_nonadditivity, steric_clearance,
)
from pincer_catmech.quantum.xtb_backend import parse_energy

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def reference_geometries():
    base = REPO / "data/reference_thermochemistry"
    return (read(base / "hemiaminal/stationarity/accepted_geometry.xyz"),
            read(base / "tert_butanol/stationarity/accepted_geometry.xyz"))


@pytest.mark.parametrize("n", [1, 2, 3])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_directional_graph_mapped_clusters(n, seed, reference_geometries):
    substrate, solvent = reference_geometries
    c = build_solvation_cluster(n, seed=seed, substrate=substrate, solvent=solvent)
    assert len(c.atoms) == 28 + 15 * n
    assert np.array_equal(c.atoms.positions[:28], substrate.positions)
    assert sorted(i for f in c.fragments for i in f) == list(range(len(c.atoms)))
    assert len(c.hydrogen_bonds) == n
    assert connectivity_check(c.atoms, c.bonds)["retained"]
    assert all(row["within_initial_window"] for row in hydrogen_bond_metrics(c.atoms, c.hydrogen_bonds))
    fragment_ids = [entry["fragment"] for entry in c.atom_map]
    assert steric_clearance(c.atoms, fragment_ids, c.hydrogen_bonds)["accepted"]
    assert c.atoms[c.substrate_sites["oxygen"]].symbol == "O"
    assert c.atoms[c.substrate_sites["nitrogen"]].symbol == "N"
    if n >= 2:
        assert {e["donor"] for e in c.hydrogen_bonds[:2]} == {c.substrate_sites["oxygen"], c.substrate_sites["nitrogen"]}
    if n == 3:
        assert c.hydrogen_bonds[2]["donor"] in c.fragments[1]
    assert not c.metadata()["is_transition_state"]


def test_seed_reproducibility_and_distinct_orientations(reference_geometries):
    a, b = reference_geometries
    one = build_solvation_cluster(2, seed=2, substrate=a, solvent=b)
    two = build_solvation_cluster(2, seed=2, substrate=a, solvent=b)
    other = build_solvation_cluster(2, seed=0, substrate=a, solvent=b)
    np.testing.assert_array_equal(one.atoms.positions, two.atoms.positions)
    assert not np.allclose(one.atoms.positions[28:], other.atoms.positions[28:])


def test_wrong_graph_identity_rejected(reference_geometries):
    a, b = reference_geometries
    bad = a.copy()
    bad.positions[0] += [10, 0, 0]
    with pytest.raises(ValueError, match="connectivity"):
        build_solvation_cluster(1, substrate=bad, solvent=b)
    bad = a.copy()
    bad.numbers[[0, 2]] = bad.numbers[[2, 0]]
    with pytest.raises(ValueError, match="atom order"):
        build_solvation_cluster(1, substrate=bad, solvent=b)


@pytest.mark.parametrize("bad", [0, 4, -1, True, 1.2])
def test_bounded_solvent_counts(bad):
    with pytest.raises(ValueError):
        build_solvation_cluster(bad)


def test_accidental_overlap_rejected():
    atoms = Atoms("OO", positions=[[0, 0, 0], [.2, 0, 0]])
    result = steric_clearance(atoms, [0, 1], [])
    assert not result["accepted"]
    assert result["clashing_pairs"] == [[0, 1]]


def test_hbond_distance_alone_does_not_pass_angle():
    # D and A lie on the same side of H: short O...H but wrong direction.
    atoms = Atoms("NHO", positions=[[0, 0, 0], [1, 0, 0], [-.85, 0, 0]])
    metric = hydrogen_bond_metrics(atoms, [{"donor": 0, "hydrogen": 1, "acceptor": 2}])[0]
    assert metric["H_acceptor_distance_A"] == pytest.approx(1.85)
    assert not metric["within_initial_window"]


def test_association_balances_all_solvent_molecules():
    assert association_energy(-42., -10., [-10.] * 3) == pytest.approx(-2.)
    assert association_energy(-42., -10., [-10.] * 2) == pytest.approx(-12.)
    with pytest.raises(ValueError):
        association_energy(float("nan"), -10., [-10.])


@pytest.mark.parametrize("n_fragments", [2, 3, 4])
def test_frozen_many_body_recovers_known_nonadditive_term(n_fragments):
    monomers = np.array([-10., -4., -5., -3.])[:n_fragments]
    pair_increments = {ij: -.1 * (1 + ij[0] + ij[1]) for ij in combinations(range(n_fragments), 2)}
    pair_energies = {ij: monomers[ij[0]] + monomers[ij[1]] + e for ij, e in pair_increments.items()}
    genuine_remainder = -.03 if n_fragments >= 3 else 0.
    total = monomers.sum() + sum(pair_increments.values()) + genuine_remainder
    result = many_body_nonadditivity(total, monomers, pair_energies)
    assert result["beyond_pair_nonadditivity_eV"] == pytest.approx(genuine_remainder, abs=1e-13)
    assert result["frozen_interaction_eV"] == pytest.approx(result["pair_sum_eV"] + genuine_remainder)
    if pair_energies:
        pair_energies.pop(next(iter(pair_energies)))
        with pytest.raises(ValueError, match="all unique"):
            many_body_nonadditivity(total, monomers, pair_energies)


def test_published_native_evidence_matches_balanced_energies():
    directory = REPO / "data/phase3/solvation"
    if not (directory / "native_evidence.zip").exists():
        pytest.skip("native physical campaign has not been executed")
    summary = json.loads((directory / "summary.json").read_text())
    manifest = json.loads((directory / "native_manifest.json").read_text())
    assert hashlib.sha256((directory / "native_evidence.zip").read_bytes()).hexdigest() == summary["native_archive"]["sha256"]
    with zipfile.ZipFile(directory / "native_evidence.zip") as archive:
        for entry in manifest:
            assert hashlib.sha256(archive.read(entry["path"])).hexdigest() == entry["sha256"]
        for row in summary["clusters"]:
            if not row["accepted_for_association"]:
                continue
            output = archive.read(f"optimizations/{row['name']}/xtb.out").decode()
            assert "GEOMETRY OPTIMIZATION CONVERGED" in output
            energy = parse_energy(output)
            assert energy == pytest.approx(row["native_result"]["energy_eV"], abs=1e-10)
            refs = summary["references"]
            expected = energy - refs["hemiaminal"]["energy_eV"] - row["n_solvent"] * refs["tert_butanol"]["energy_eV"]
            assert expected == pytest.approx(row["association_energy_eV"], abs=1e-10)
            mapping = json.loads((directory / "clusters" / row["name"] / "mapping.json").read_text())
            assert all(x["within_initial_window"] for x in mapping["initial_hbond_metrics"])
    for selected in summary["selected"]:
        stationary = json.loads((directory / "clusters" / selected["name"] / "stationarity/result.json").read_text())
        if stationary["accepted"]:
            attempt = stationary["attempts"][stationary["accepted_attempt"]]
            assert attempt["imaginary_count"] == 0
            assert attempt["external_modes_projected"] == 6
