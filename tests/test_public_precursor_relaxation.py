"""Offline scientific acceptance and immutable-output boundaries; no QC execution."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase.io import read


@pytest.fixture
def cli():
    spec = importlib.util.spec_from_file_location("public_relax", Path(__file__).resolve().parents[1] / "scripts/run_public_precursor_relaxation.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def seed(cli):
    p = cli.SOURCE / "Mn1"
    return read(p / "crystal_candidate.xyz"), json.loads((p / "metadata.json").read_text())


def test_original_source_graph_accepted(cli, seed):
    atoms, metadata = seed
    result = cli.graph_identity(atoms, atoms, metadata)
    assert result["retained"] and result["essential_bonds_retained"]
    assert result["NH_distance_angstrom"] < .9


@pytest.mark.parametrize("site", ["halide", "CO", "NH", "P"])
def test_loss_of_named_precursor_component_rejected(cli, seed, site):
    atoms, metadata = seed
    moved = atoms.copy()
    index = {"halide": metadata["halide_indices"][0], "CO": metadata["carbonyl_indices"][0][1],
             "NH": metadata["ligand_proton_index"], "P": metadata["donor_indices"][-1]}[site]
    moved.positions[index] += [20, 0, 0]
    assert cli.graph_identity(moved, atoms, metadata)["retained"] is False


def test_atom_order_mismatch_rejected(cli, seed):
    atoms, metadata = seed
    moved = atoms.copy()
    moved.numbers[[0, 1]] = moved.numbers[[1, 0]]
    assert not cli.graph_identity(moved, atoms, metadata)["retained"]


@pytest.mark.parametrize("minimum,asymmetry,accepted", [(1., .001, True), (-.01, .001, False), (0., .001, False), (1., .1, False)])
def test_full_hessian_acceptance_keeps_all_negative_and_zero_modes(cli, minimum, asymmetry, accepted):
    spectrum = SimpleNamespace(frequencies_cm1=np.array([minimum] + [10.] * 179), antisymmetry_relative=asymmetry)
    assert cli.validate_spectrum(spectrum, 62)["hessian_minimum_pass"] is accepted


def test_partial_hessian_rejected(cli):
    with pytest.raises(ValueError, match="cardinality"):
        cli.validate_spectrum(SimpleNamespace(frequencies_cm1=np.ones(30), antisymmetry_relative=0.), 62)


def test_existing_and_nested_paths_rejected(cli, tmp_path):
    existing = tmp_path / "old"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        cli.fresh_paths(existing, tmp_path / "scratch")
    with pytest.raises(ValueError):
        cli.fresh_paths(tmp_path / "fresh", tmp_path / "fresh/sub")


def test_total_and_individual_budgets_both_apply(cli, monkeypatch):
    monkeypatch.setattr(cli.time, "time", lambda: 100.)
    assert cli.remaining_budget({"quantum_wall_seconds": 400.}, 1000.) == 20.
    assert cli.remaining_budget({"quantum_wall_seconds": 0.}, 105.) == 5.
    assert cli.remaining_budget({"quantum_wall_seconds": 0.}, 90.) == 0.


def test_preserve_all_native_files(cli, tmp_path):
    work, folder = tmp_path / "scratch", tmp_path / "out"
    work.mkdir()
    folder.mkdir()
    (work / "xtbrestart").write_bytes(b"MOCK RESTART")
    (work / "xtb.out").write_text("MOCK; NO REAL CALCULATION")
    receipt = cli.preserve_native(work, folder)
    assert receipt["file_count"] == 2
    assert (work / "xtbrestart").exists()
    assert cli.digest(folder / "native/xtbrestart") == cli.digest(work / "xtbrestart")
