"""Public-CIF provenance, crystallographic identity and failure-boundary tests."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms


@pytest.fixture(scope="module")
def cli():
    spec = importlib.util.spec_from_file_location("public_catalyst", Path(__file__).resolve().parents[1] / "scripts/import_public_catalyst.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def original():
    return Path(__file__).resolve().parents[1] / "data/phase4/public_structure/original.cif"


def test_preamble_repair_changes_no_scientific_fields(cli, original):
    raw = original.read_text()
    normalized, receipt = cli.normalized_cif(raw)
    assert receipt["changed"]
    assert normalized.split("\n", 1)[1] == raw[raw.index("_audit_creation_method"):]
    assert cli.digest(original) == cli.EXPECTED_SHA256


def test_unknown_preamble_fails_closed(cli):
    with pytest.raises(ValueError, match="Unsupported"):
        cli.normalized_cif("Some other structure\n_audit_creation_method 'x'")


def test_existing_data_header_left_unchanged(cli):
    data = "data_x\n_tag value\n"
    assert cli.normalized_cif(data) == (data, {"changed": False})


def test_formula_unknown_syntax_rejected(cli):
    with pytest.raises(ValueError):
        cli.formula_counts("'C18 H37 Br Mn N O2 P2 + solvent'")


def test_real_cif_occupancy_hydrogens_and_bonds(cli, original):
    parsed = cli.parse_structure(cli.normalized_cif(original.read_text())[0])
    assert len(parsed["sites"]) == 124 and parsed["Z"] == 8
    assert parsed["formula"] == cli.EXPECTED_FORMULA
    assert all(s["occupancy"] == 1 and s["disorder_group"] == "." for s in parsed["sites"])
    assert sum(s["symbol"] == "H" for s in parsed["sites"]) == 74
    errors = cli.check_published_bonds(parsed)
    assert max(x["absolute_deviation_angstrom"] for x in errors) < .002


def test_periodic_molecule_unwrapped_without_duplication(cli):
    atoms = Atoms("HH", scaled_positions=[[.97, .5, .5], [.04, .5, .5]], cell=[10, 10, 10], pbc=True)
    components, graph, infinite = cli.periodic_components(atoms)
    assert len(components) == 1 and not infinite
    c = components[0]
    unwrapped = atoms.positions[c["indices"]] + np.array(c["offsets"]) @ atoms.cell.array
    assert np.linalg.norm(unwrapped[1] - unwrapped[0]) == pytest.approx(.7)


def test_infinite_periodic_chain_detected(cli):
    atoms = Atoms("C", scaled_positions=[[0, 0, 0]], cell=[1.4, 10, 10], pbc=True)
    components, graph, infinite = cli.periodic_components(atoms)
    assert infinite


def test_substituted_download_rejected_before_output(cli, tmp_path):
    source = tmp_path / "fake.cif"
    source.write_text("bad")
    with pytest.raises(ValueError, match="checksum"):
        cli.import_structure(source, tmp_path / "result")
    assert not (tmp_path / "result").exists()


def test_output_never_overwrites(cli, original, tmp_path):
    out = tmp_path / "existing"
    out.mkdir()
    (out / "prior").write_text("keep")
    with pytest.raises(FileExistsError):
        cli.import_structure(original, out)
    assert (out / "prior").read_text() == "keep"


def test_complete_real_import_is_candidate_not_solution_or_spin_proof(cli, original, tmp_path):
    out = tmp_path / "fresh_import"
    report = cli.import_structure(original, out)
    assert report["status"] == "structure_extraction_verified"
    assert report["unit_cell_atom_count"] == 496
    assert report["periodic_component_count"] == 8
    assert len(report["molecules"]) == 2
    assert report["quantum_submission_ready"] is False
    assert report["charge_spin_experimentally_verified"] is False
    metadata = json.loads((out / "Mn1/metadata.json").read_text())
    assert metadata["charge"] == 0 and metadata["multiplicity"] == 1
    assert metadata["proton_site_element"] == "N"
    assert len(metadata["carbonyl_indices"]) == 2 and len(metadata["halide_indices"]) == 1
    assert metadata["structure_origin"] == "public_crystal_derived_candidate"
    assert all(metadata["provenance"][key] for key in ("source_url", "source_identifier", "structure_derivation"))
    assert metadata["third_party_data_license"]["name"] == "CC BY-NC 4.0"
    assert report["molecules"][0]["graph_matches_published_bonds"]
    manifest = json.loads((out / "file_manifest.json").read_text())
    assert all(cli.digest(out / path) == checksum for path, checksum in manifest.items())
