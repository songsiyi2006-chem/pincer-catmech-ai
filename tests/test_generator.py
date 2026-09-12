"""Chemical graph, atom-mapping and coordinate checks, never energy mocks."""

from collections import Counter
import json

import numpy as np
import pytest

from pincer_catmech.generators import CatalystSpec, enumerate_specs, generate_catalyst, write_library


def test_matrix_is_exactly_twenty_four_distinct_designs():
    specs = enumerate_specs()
    assert len(specs) == len({s.catalyst_id for s in specs}) == 24
    assert Counter(s.metal for s in specs) == {"Ru": 6, "Mn": 6, "Fe": 6, "Co": 6}
    assert Counter(s.substituent for s in specs) == {"Ph": 12, "iPr": 12}


@pytest.fixture(scope="module")
def structures():
    pytest.importorskip("rdkit", reason="RDKit is required for actual ligand graph/3D generation")
    return {
        (spec.catalyst_id, state): generate_catalyst(spec, state, seed=20260913 + index)
        for index, spec in enumerate(enumerate_specs())
        for state in ("active", "hydrogenated", "protonated_reference")
    }


def test_all_seventy_two_have_finite_nonoverlapping_coordinates(structures):
    assert len(structures) == 72
    for generated in structures.values():
        assert 30 < len(generated.atoms) < 100
        assert np.isfinite(generated.atoms.positions).all()
        distances = generated.atoms.get_all_distances() + np.eye(len(generated.atoms)) * 999
        assert distances.min() > .75
        assert all(1.8 < generated.atoms.get_distance(0, i) < 2.7 for i in generated.donor_indices)


def test_state_transformations_conserve_nonhydrogen_atom_inventory(structures):
    for spec in enumerate_specs():
        active, hydrogenated, protonated = [structures[(spec.catalyst_id, s)] for s in ("active", "hydrogenated", "protonated_reference")]
        counts = [Counter(g.atoms.get_chemical_symbols()) for g in (active, hydrogenated, protonated)]
        assert counts[1] - counts[0] == {"H": 2}
        assert counts[2] - counts[0] == {"H": 1}
        assert not (counts[0] - counts[1])
        assert not (counts[0] - counts[2])
        assert (active.charge, hydrogenated.charge, protonated.charge) == (0, 0, 1)
        assert active.proton_site_index == hydrogenated.proton_site_index == protonated.proton_site_index
        assert active.donor_indices == hydrogenated.donor_indices == protonated.donor_indices
        assert active.ligand_proton_index is None
        assert hydrogenated.ligand_proton_index is not None
        assert len(hydrogenated.hydride_indices) == len(active.hydride_indices) + 1


def test_formal_graph_charges_and_electron_occupation_parity(structures):
    from rdkit import Chem
    for generated in structures.values():
        ligand = Chem.MolFromSmiles(generated.ligand_smiles)
        assert ligand is not None
        assert Chem.GetFormalCharge(ligand) == (-1 if generated.state == "active" else 0)
        assert sum(atom.GetNumRadicalElectrons() for atom in ligand.GetAtoms()) == 0
        assert (sum(generated.atoms.numbers) - generated.charge - generated.unpaired_electrons) % 2 == 0
        assert generated.metadata()["spin_provenance"].startswith("Nominal")


def test_backbones_have_explicit_n_c_o_proton_sites_and_correct_donors(structures):
    for generated in structures.values():
        symbols = generated.atoms.get_chemical_symbols()
        expected = {"macho_pnp": ("N", ["P", "N", "P"]),
                    "pyridine_pnn": ("C", ["P", "N", "N"]),
                    "bipyridine_pnnoh": ("O", ["P", "N", "N"])}[generated.spec.backbone]
        assert symbols[generated.proton_site_index] == expected[0]
        assert [symbols[i] for i in generated.donor_indices] == expected[1]
        assert generated.metadata()["P_M_P_applicable"] == (generated.spec.backbone == "macho_pnp")
        if generated.spec.backbone == "bipyridine_pnnoh":
            assert 4 in generated.source_atom_maps and 5 in generated.source_atom_maps
            assert symbols[generated.source_atom_maps[4]] == "C"
            assert symbols[generated.source_atom_maps[5]] == "O"


def test_ancillary_sets_are_explicit_and_no_energies_fabricated(structures):
    for generated in structures.values():
        assert len(generated.carbonyl_indices) == (2 if generated.spec.metal == "Mn" else 1)
        baseline_hydrides = 1 if generated.spec.metal in ("Ru", "Fe") else 0
        assert len(generated.hydride_indices) == baseline_hydrides + (generated.state == "hydrogenated")
        for c, o in generated.carbonyl_indices:
            assert generated.atoms[c].symbol == "C" and generated.atoms[o].symbol == "O"
            assert generated.atoms.get_distance(c, o) == pytest.approx(1.14)
        metadata = generated.metadata()
        assert not any(key in metadata for key in ("energy", "electronic_energy", "barrier", "global_minimum"))
        assert metadata["structure_status"].startswith("Designed starting geometry")
        assert json.loads(json.dumps(metadata))["catalyst_id"] == generated.spec.catalyst_id


def test_invalid_identity_and_state_are_rejected():
    with pytest.raises(ValueError):
        CatalystSpec("Pt", "macho_pnp", "Ph")
    pytest.importorskip("rdkit")
    with pytest.raises(ValueError, match="State"):
        generate_catalyst(enumerate_specs()[0], state="unknown")


def test_seeded_coordinates_are_reproducible():
    pytest.importorskip("rdkit")
    spec = CatalystSpec("Mn", "pyridine_pnn", "iPr")
    first = generate_catalyst(spec, seed=718)
    second = generate_catalyst(spec, seed=718)
    np.testing.assert_allclose(first.atoms.positions, second.atoms.positions, atol=1e-10)


def test_torsion_seed_changes_substituent_without_moving_donors():
    pytest.importorskip("rdkit")
    spec = CatalystSpec("Ru", "macho_pnp", "iPr")
    first = generate_catalyst(spec, seed=71)
    second = generate_catalyst(spec, seed=71, arm_torsion_deg=30)
    np.testing.assert_allclose(first.atoms.positions[list(first.donor_indices)], second.atoms.positions[list(second.donor_indices)], atol=1e-10)
    assert not np.allclose(first.atoms.positions, second.atoms.positions)
