"""Provenance, fixed-temperature selection and stoichiometric software tests."""
import hashlib
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

from ase import Atoms
import pytest


@pytest.fixture
def exporter():
    path = Path(__file__).resolve().parents[1] / "scripts/export_thermochemistry_dataset.py"
    specification = importlib.util.spec_from_file_location("thermochemistry_export_under_test", path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def test_clone_prefers_local_same_hash_and_rejects_tampering(tmp_path, monkeypatch, exporter):
    monkeypatch.setattr(exporter, "REPO", tmp_path)
    folder = tmp_path / "certified"
    folder.mkdir()
    artifact = folder / "hessian_00.npz"
    # Bytes test only: this deliberately is not a numerical research Hessian.
    artifact.write_bytes(b"software integrity fixture")
    expected = hashlib.sha256(artifact.read_bytes()).hexdigest()
    resolved = exporter.resolve_evidence(folder / "result.json", r"C:\unavailable\old-checkout\hessian_00.npz", expected)
    assert resolved == artifact
    artifact.write_bytes(b"changed")
    with pytest.raises(ValueError, match="matching SHA256"):
        exporter.resolve_evidence(folder / "result.json", r"C:\unavailable\old-checkout\hessian_00.npz", expected)


def test_one_certificate_selected_at_base_temperature_for_entire_grid(tmp_path, monkeypatch, exporter):
    monkeypatch.setattr(exporter, "REPO", tmp_path)
    # Deliberate abstract comparison values; no model or molecule is simulated.
    first = SimpleNamespace(key="same_species", certificate=tmp_path / "a.json",
                            rows={383.15: {"G_298_qRRHO_sol": 1}, 440.: {"G_298_qRRHO_sol": 3}})
    second = SimpleNamespace(key="same_species", certificate=tmp_path / "b.json",
                             rows={383.15: {"G_298_qRRHO_sol": 2}, 440.: {"G_298_qRRHO_sol": 0}})
    selected = exporter.select_certificates([second, first])
    assert selected["same_species"] is first
    assert selected["same_species"].rows[440.]["G_298_qRRHO_sol"] == 3


def test_reaction_balance_checks_charge_as_well_as_every_element(exporter):
    proton = SimpleNamespace(atoms=Atoms("H"), charge=1)
    hydroxide = SimpleNamespace(atoms=Atoms("OH"), charge=-1)
    water = SimpleNamespace(atoms=Atoms("H2O"), charge=0)
    assert exporter.reaction_balance([proton, hydroxide], [water])["balanced"]
    water.charge = 1
    assert not exporter.reaction_balance([proton, hydroxide], [water])["balanced"]
    assert not exporter.reaction_balance([proton], [hydroxide])["balanced"]


def test_malformed_certificate_preserved_as_rejected_ledger_entry(tmp_path, monkeypatch, exporter):
    monkeypatch.setattr(exporter, "REPO", tmp_path)
    folder = tmp_path / "catalyst__active"
    folder.mkdir()
    certificate = folder / "result.json"
    certificate.write_text('{"incomplete":', encoding="utf-8")
    ledger, accepted = exporter.validate_certificate(certificate, {})
    assert accepted is None and not ledger["eligible"]
    assert ledger["key"] == "catalyst|active"
    assert ledger["validation_status"] == "malformed_certificate"
    assert ledger["certificate_sha256"] == hashlib.sha256(certificate.read_bytes()).hexdigest()


def test_cycle_species_cancellation_does_not_confuse_equal_formulas(exporter):
    # Abstract stoichiometry only; no simulated or research energies.
    steps = [(["cat", "alcohol"], ["catH2", "aldehyde"]),
             (["aldehyde", "aniline"], ["hemiaminal"]),
             (["hemiaminal"], ["imine", "water"]),
             (["catH2", "imine"], ["cat", "product_amine"])]
    net = [(["alcohol", "aniline"], ["product_amine", "water"])]
    assert exporter.species_stoichiometry(steps) == exporter.species_stoichiometry(net)
    # Consuming another isomer of the same formula must not cancel this imine.
    steps[-1] = (["catH2", "other_imine_isomer"], ["cat", "product_amine"])
    assert exporter.species_stoichiometry(steps) != exporter.species_stoichiometry(net)
