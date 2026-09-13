"""Atom-conservation and evidence-boundary checks; no synthetic QC results."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("phase3_wire_driver", REPO / "scripts" / "run_phase3_proton_wire.py")
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


@pytest.mark.parametrize("seed,water_distance", [(2, 4.3), (1, 5.4)])
def test_actual_source_conserves_all_atom_identities_and_transfers_distinct_protons(seed, water_distance):
    source = REPO / "data" / "phase3" / "solvation" / "clusters" / f"hemiaminal_tbuoh_1_seed_{seed:02d}"
    reactant, product, graphs, coordinates, cn, metadata = DRIVER.mapped_endpoint_seed(source, water_distance)
    assert len(reactant) == len(product) == 43
    np.testing.assert_array_equal(reactant.numbers, product.numbers)
    assert metadata["reactant_seed_identity"]["retained"]
    sites = metadata["mapped_indices"]
    def neighbors(graph, atom):
        return {j if i == atom else i for i, j in graph if atom in (i, j)}
    assert neighbors(graphs[0], sites["N_H"]) == {sites["N"]}
    assert neighbors(graphs[0], sites["shuttle_H"]) == {sites["shuttle_O"]}
    assert neighbors(graphs[1], sites["N_H"]) == {sites["shuttle_O"]}
    assert neighbors(graphs[1], sites["shuttle_H"]) == {sites["leaving_O"]}
    assert neighbors(graphs[1], sites["leaving_O"]) == {sites["original_O_H"], sites["shuttle_H"]}
    assert 1.15 <= product.get_distance(*cn) <= 1.42
    assert product.get_distance(sites["C"], sites["leaving_O"]) == pytest.approx(water_distance)
    assert len(coordinates) == 4
    assert metadata["closed_six_or_eight_membered_TS_demonstrated"] is False
    assert metadata["seed_geometry_is_stationary_point_evidence"] is False


def test_rejects_ionic_model_before_using_neutral_wire(tmp_path):
    source = REPO / "data" / "phase3" / "solvation" / "clusters" / "hemiaminal_tbuoh_1_seed_02" / "mapping.json"
    metadata = json.loads(source.read_text(encoding="utf-8"))
    metadata["charge"] = -1
    (tmp_path / "mapping.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="neutral"):
        DRIVER.mapped_endpoint_seed(tmp_path)


def test_rejects_corrupted_map_before_preparing_proton_transfers(tmp_path):
    source = REPO / "data" / "phase3" / "solvation" / "clusters" / "hemiaminal_tbuoh_1_seed_02"
    metadata = json.loads((source / "mapping.json").read_text(encoding="utf-8"))
    metadata["atom_map"][17]["global_index"] = 42
    (tmp_path / "mapping.json").write_text(json.dumps(metadata), encoding="utf-8")
    (tmp_path / "optimized.xyz").write_bytes((source / "optimized.xyz").read_bytes())
    with pytest.raises(ValueError, match="complete ordered atom mapping"):
        DRIVER.mapped_endpoint_seed(tmp_path)
