"""Atom/bond/geometry checks for endpoint construction, without mock energies."""

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from pincer_catmech.generators import CatalystSpec, generate_catalyst

SPEC = importlib.util.spec_from_file_location("pincer_neb_driver", Path(__file__).resolve().parents[1] / "scripts/run_neb_campaign.py")
DRIVER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DRIVER
SPEC.loader.exec_module(DRIVER)


@pytest.fixture(scope="module")
def assembly():
    pytest.importorskip("rdkit")
    catalyst = generate_catalyst(CatalystSpec("Ru", "macho_pnp", "iPr"), seed=122)
    return DRIVER.assemble_endpoints(catalyst.atoms, catalyst.metadata(), seed=641)


def test_complete_atom_mapping_and_electronic_state_are_conserved(assembly):
    np.testing.assert_array_equal(assembly.reactant.numbers, assembly.product.numbers)
    np.testing.assert_array_equal(assembly.reactant.get_masses(), assembly.product.get_masses())
    assert assembly.reactant.get_chemical_formula() == assembly.product.get_chemical_formula()
    for key in ("charge", "unpaired", "uhf"):
        assert assembly.reactant.info[key] == assembly.product.info[key] == 0
    assert assembly.metadata["complete_atom_count"] == assembly.metadata["catalyst_atom_count"] + 16


def test_hydrogens_remain_identified_and_bracket_both_transfers(assembly):
    transfer = assembly.transfer_indices
    assert transfer["proton"] != transfer["hydride"]
    for channel in ("proton", "hydride"):
        hydrogen, donor, acceptor = (transfer[channel], transfer[channel + "_donor"], transfer[channel + "_acceptor"])
        assert assembly.reactant[hydrogen].symbol == "H"
        reactant_progress = assembly.reactant.get_distance(hydrogen, donor) - assembly.reactant.get_distance(hydrogen, acceptor)
        product_progress = assembly.product.get_distance(hydrogen, donor) - assembly.product.get_distance(hydrogen, acceptor)
        assert reactant_progress < -.1
        assert product_progress > .30
    assert assembly.reactant[transfer["proton_donor"]].symbol == "O"
    assert assembly.reactant[transfer["hydride_donor"]].symbol == "C"


def test_alcohol_to_aldehyde_carbonyl_geometry_without_clashes(assembly):
    o = assembly.metadata["substrate_oxygen_index"]
    c = assembly.metadata["substrate_carbonyl_carbon_index"]
    assert 1.35 < assembly.reactant.get_distance(c, o) < 1.55
    assert 1.15 < assembly.product.get_distance(c, o) < 1.30
    for endpoint in (assembly.reactant, assembly.product):
        assert np.min(endpoint.get_all_distances() + np.eye(len(endpoint)) * 999) > .70
    assert "NOT an energy" in assembly.metadata["placement_score_units"]


def test_protonated_reference_cannot_masquerade_as_active_endpoint():
    pytest.importorskip("rdkit")
    catalyst = generate_catalyst(CatalystSpec("Mn", "macho_pnp", "iPr"), state="protonated_reference")
    with pytest.raises(ValueError, match="deprotonated active"):
        DRIVER.assemble_endpoints(catalyst.atoms, catalyst.metadata())


def test_ligand_carbon_hydrogen_migration_is_an_identity_failure():
    pytest.importorskip("rdkit")
    catalyst = generate_catalyst(CatalystSpec("Ru", "macho_pnp", "iPr"), seed=122)
    assert DRIVER.ligand_identity_screen(catalyst.atoms, catalyst.metadata())["ligand_retained"]
    carbon, hydrogen = next((i, j) if catalyst.atoms[i].symbol == "C" else (j, i)
        for i, j in catalyst.ligand_bond_indices
        if {catalyst.atoms[i].symbol, catalyst.atoms[j].symbol} == {"C", "H"})
    altered = catalyst.atoms.copy()
    altered.positions[hydrogen] = altered.positions[0] + np.array([0., 0., -1.55])
    screen = DRIVER.ligand_identity_screen(altered, catalyst.metadata())
    assert not screen["retained"]
    assert not screen["ligand_retained"]
    assert tuple(sorted((carbon, hydrogen))) in screen["lost_bonds"]


def test_restart_mode_selection_uses_full_matrix_and_binds_exact_geometry(tmp_path):
    """Manufactured curvature tests the adapter; it is not quantum campaign data."""
    from ase import Atoms
    atoms = Atoms("OHNCHRu", positions=[[0, 0, 0], [1, 0, 0], [2, .3, 0],
                  [0, 2, 0], [1, 2, .2], [2, 2, 0]])
    transfer = dict(proton=1, proton_donor=0, proton_acceptor=2,
                    hydride=4, hydride_donor=3, hydride_acceptor=5)
    source = tmp_path / "manufactured_hessian.npz"
    def save(sign):
        np.savez(source, hessian_ev_a2=sign * np.eye(3 * len(atoms)),
                 positions=atoms.positions, atomic_numbers=atoms.numbers,
                 masses=atoms.get_masses(), displacement_angstrom=.005)
    save(-1)
    mode, provenance = DRIVER.load_restart_hessian_mode(atoms, source, transfer)
    assert mode.shape == atoms.positions.shape
    assert np.isfinite(mode).all() and np.linalg.norm(mode) > 0
    assert provenance["selected_mode"]["combined_overlap"] == max(
        row["combined_overlap"] for row in provenance["negative_modes"])
    assert len(provenance["source_sha256"]) == 64
    shifted = atoms.copy()
    shifted.positions[0, 0] += .001
    with pytest.raises(ValueError, match="exact candidate"):
        DRIVER.load_restart_hessian_mode(shifted, source, transfer)
    save(1)
    with pytest.raises(ValueError, match="no actual negative"):
        DRIVER.load_restart_hessian_mode(atoms, source, transfer)
