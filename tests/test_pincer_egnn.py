"""Manufactured tensors verify software symmetries; none are chemical labels."""
import math

import pytest
import torch

from pincer_catmech.models.pincer_egnn import (
    PincerEGNN, MaskedTargetNormalizer, complete_edge_index, grouped_split, masked_mse_loss,
)


@pytest.fixture
def sample():
    torch.manual_seed(712)
    model = PincerEGNN(5, graph_feature_dim=2, hidden_dim=16, num_layers=3).double()
    return model, torch.randn(9, 5, dtype=torch.float64), torch.randn(9, 3, dtype=torch.float64)


def euler_rotation():
    a, b, c = .731, -1.247, 2.091
    cx, sx, cy, sy, cz, sz = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    rx = torch.tensor([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], dtype=torch.float64)
    ry = torch.tensor([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], dtype=torch.float64)
    rz = torch.tensor([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], dtype=torch.float64)
    return rz @ ry @ rx


@pytest.mark.parametrize("parity", [1., -1.])
def test_so3_scalar_and_se3_coordinate_equivariance_with_o3_parity(sample, parity):
    model, features, positions = sample
    rotation = parity * euler_rotation()
    torch.testing.assert_close(rotation.T @ rotation, torch.eye(3, dtype=torch.float64), atol=1e-14, rtol=0)
    context = torch.tensor([[.37, -1.2]], dtype=torch.float64)
    original = model(features, positions, graph_features=context)
    translation = torch.tensor([5.37, -2.19, 8.21], dtype=torch.float64)
    moved = model(features, positions @ rotation.T + translation, graph_features=context)
    assert (original.positions - positions).abs().max() > 1e-8  # Coordinate update actually operates.
    assert (original.scalars - moved.scalars).abs().max() < 1e-6  # Requested eV acceptance threshold.
    torch.testing.assert_close(original.scalars, moved.scalars, atol=2e-12, rtol=0)
    torch.testing.assert_close(original.positions @ rotation.T + translation, moved.positions, atol=2e-12, rtol=0)
    torch.testing.assert_close(original.node_features, moved.node_features, atol=2e-12, rtol=0)


def test_permutation_of_nodes_edges_and_edge_features():
    torch.manual_seed(401)
    model = PincerEGNN(4, edge_feature_dim=2, hidden_dim=12, num_layers=2).double()
    features, positions = torch.randn(7, 4, dtype=torch.float64), torch.randn(7, 3, dtype=torch.float64)
    edges = complete_edge_index(torch.zeros(7, dtype=torch.long))
    attributes = torch.randn(edges.shape[1], 2, dtype=torch.float64)
    original = model(features, positions, edge_index=edges, edge_features=attributes)
    permutation = torch.tensor([5, 1, 0, 6, 3, 2, 4])
    inverse = torch.argsort(permutation)
    edge_permutation = torch.randperm(edges.shape[1])
    reordered = model(features[permutation], positions[permutation],
                      edge_index=inverse[edges[:, edge_permutation]], edge_features=attributes[edge_permutation])
    torch.testing.assert_close(original.scalars, reordered.scalars, atol=2e-12, rtol=0)
    torch.testing.assert_close(original.positions[permutation], reordered.positions, atol=2e-12, rtol=0)
    torch.testing.assert_close(original.node_features[permutation], reordered.node_features, atol=2e-12, rtol=0)


def test_batch_independence_including_empty_graphs_and_graph_permutation(sample):
    model, features, positions = sample
    batch = torch.tensor([0] * 4 + [2] * 5)
    context = torch.randn(4, 2, dtype=torch.float64)
    packed = model(features, positions, batch=batch, num_graphs=4, graph_features=context)
    assert packed.graph_mask.tolist() == [True, False, True, False]
    assert torch.equal(packed.scalars[[1, 3]], torch.zeros(2, 3, dtype=torch.float64))
    for graph, indices in ((0, slice(0, 4)), (2, slice(4, 9))):
        single = model(features[indices], positions[indices], graph_features=context[graph:graph + 1])
        torch.testing.assert_close(single.scalars[0], packed.scalars[graph], atol=2e-12, rtol=0)
        torch.testing.assert_close(single.positions, packed.positions[indices], atol=2e-12, rtol=0)
    graph_order = torch.tensor([2, 3, 0, 1])
    permuted = model(features, positions, batch=graph_order[batch], num_graphs=4, graph_features=context[graph_order])
    torch.testing.assert_close(permuted.scalars[graph_order], packed.scalars, atol=2e-12, rtol=0)


def test_padding_mask_cannot_change_active_predictions(sample):
    model, features, positions = sample
    context = torch.zeros(1, 2, dtype=torch.float64)
    reference = model(features, positions, graph_features=context)
    padded_features = torch.cat((features, torch.full((3, 5), 999., dtype=torch.float64)))
    padded_positions = torch.cat((positions, torch.full((3, 3), -321., dtype=torch.float64)))
    output = model(padded_features, padded_positions, node_mask=torch.tensor([True] * 9 + [False] * 3), graph_features=context)
    torch.testing.assert_close(reference.scalars, output.scalars, atol=2e-12, rtol=0)
    torch.testing.assert_close(reference.positions, output.positions[:9], atol=2e-12, rtol=0)
    assert torch.equal(output.positions[9:], padded_positions[9:])
    assert torch.count_nonzero(output.node_features[9:]) == 0


@pytest.mark.parametrize("nodes,graphs", [(0, 0), (0, 3), (1, 1), (5, 2)])
def test_empty_graphs_zero_edges_and_coincident_coordinates_are_finite(nodes, graphs):
    model = PincerEGNN(3, hidden_dim=8, num_layers=2).double()
    features = torch.randn(nodes, 3, dtype=torch.float64, requires_grad=True)
    positions = torch.zeros(nodes, 3, dtype=torch.float64, requires_grad=True)
    output = model(features, positions, edge_index=torch.empty(2, 0, dtype=torch.long), num_graphs=graphs)
    assert torch.isfinite(output.scalars).all() and torch.isfinite(output.positions).all()
    torch.testing.assert_close(output.positions, positions)
    if nodes:
        interacting = model(features, positions, num_graphs=graphs)
        loss = interacting.scalars.square().sum() + interacting.positions.square().sum()
        loss.backward()
        assert torch.isfinite(features.grad).all() and torch.isfinite(positions.grad).all()
        assert all(parameter.grad is None or torch.isfinite(parameter.grad).all() for parameter in model.parameters())


def test_scalar_coordinate_gradients_are_equivariant(sample):
    model, features, positions = sample
    positions.requires_grad_()
    rotation = euler_rotation()
    moved = (positions.detach() @ rotation.T + 3.12).requires_grad_()
    first = model(features, positions).scalars[:, 0].sum()
    second = model(features, moved).scalars[:, 0].sum()
    gradient = torch.autograd.grad(first, positions)[0]
    moved_gradient = torch.autograd.grad(second, moved)[0]
    torch.testing.assert_close(gradient @ rotation.T, moved_gradient, atol=2e-12, rtol=0)


def test_missing_labels_are_masked_without_nan_or_fabricated_zero_targets():
    predicted = torch.tensor([[1., 2., 3.]], dtype=torch.float64, requires_grad=True)
    target = torch.tensor([[2., float("nan"), float("nan")]], dtype=torch.float64)
    mask = torch.tensor([[True, False, False]])
    loss = masked_mse_loss(predicted, target, mask)
    assert loss.item() == 1
    loss.backward()
    assert predicted.grad.tolist() == [[-2., 0., 0.]]
    nan_predictions = torch.full((1, 3), float("nan"), requires_grad=True)
    empty_loss = masked_mse_loss(nan_predictions, target.float(), torch.zeros_like(mask))
    assert empty_loss.item() == 0
    empty_loss.backward()
    assert torch.equal(nan_predictions.grad, torch.zeros_like(nan_predictions))
    with pytest.raises(ValueError, match="finite"):
        masked_mse_loss(predicted, target, torch.ones_like(mask))


def test_training_only_normalization_and_family_split_prevent_conformer_leakage():
    groups = [f"ligand_{group}" for group in range(6) for _ in range(4)]
    split = grouped_split(groups)
    families = {name: {groups[i] for i in indices} for name, indices in split.items()}
    assert not families["train"] & families["test"]
    assert not families["train"] & families["validation"]
    assert not families["test"] & families["validation"]
    values = torch.arange(24, dtype=torch.float64)[:, None].expand(24, 3).clone()
    values[:, 1] = float("nan")
    mask = torch.isfinite(values)
    fitted = MaskedTargetNormalizer.fit(values, mask, split["train"])
    altered = values.clone()
    altered[split["test"] + split["validation"], 0] = 1e12
    repeated = MaskedTargetNormalizer.fit(altered, mask, split["train"])
    torch.testing.assert_close(fitted.mean, repeated.mean)
    torch.testing.assert_close(fitted.scale, repeated.scale)
    assert fitted.supported.tolist() == [True, False, True]
    with pytest.raises(ValueError, match="three"):
        grouped_split(["same_family"] * 100)


def test_untrained_heads_and_invalid_edges_are_rejected(sample):
    model, features, positions = sample
    with pytest.raises(ValueError, match="untrained"):
        model.predict_trained("activation_barrier_ev", features, positions)
    with pytest.raises(ValueError, match="Cross-graph"):
        model(features, positions, batch=torch.arange(9), edge_index=torch.tensor([[0], [1]]))
    model.mark_heads_trained(["aux_conformer_relative_electronic_energy_ev"])
    assert model.trained_heads.tolist() == [False, False, True]
    assert torch.isfinite(model.predict_trained("aux_conformer_relative_electronic_energy_ev", features, positions)).all()


def test_pair_only_training_cannot_be_misrepresented_as_absolute_energy(sample):
    model, features, positions = sample
    name = "aux_conformer_relative_electronic_energy_ev"
    model.mark_heads_trained([name], pair_only=True)
    first = {"node_features": features, "positions": positions}
    second = {"node_features": features, "positions": positions * 1.03}
    with pytest.raises(ValueError, match="pair"):
        model.predict_trained(name, features, positions)
    forward = model.predict_pair_trained(name, first, second)
    reverse = model.predict_pair_trained(name, second, first)
    torch.testing.assert_close(forward, -reverse, atol=0, rtol=0)
    assert torch.count_nonzero(model.predict_pair_trained(name, first, first)) == 0


def test_source_hash_verification_accepts_only_documented_windows_newline_change():
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location("train_phase3_egnn_test", Path(__file__).resolve().parents[1] / "scripts/train_phase3_egnn.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    geometry = b"1\nreal source\nH 0 0 0\n"
    output_lf = b"TOTAL ENERGY -1.000 Eh\nnormal termination of xtb\n"
    module.verify_native_hashes(geometry, output_lf.replace(b"\n", b"\r\n"), module.digest(geometry), module.digest(output_lf))
    with pytest.raises(ValueError, match="stdout"):
        module.verify_native_hashes(geometry, output_lf.replace(b"-1.000", b"-2.000"), module.digest(geometry), module.digest(output_lf))
    with pytest.raises(ValueError, match="Geometry"):
        module.verify_native_hashes(geometry + b" ", output_lf, module.digest(geometry), module.digest(output_lf))


def test_empty_dataset_produces_explicit_untrained_readiness(tmp_path):
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location("train_phase3_egnn_empty_test", Path(__file__).resolve().parents[1] / "scripts/train_phase3_egnn.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.write_jsonl(tmp_path / "structures.jsonl", [])
    module.write_jsonl(tmp_path / "pairs.jsonl", [])
    assert (tmp_path / "pairs.jsonl").read_bytes() == b""
    result = module.run_training(tmp_path)
    assert result["status"] == "insufficient_auxiliary_data"
    assert not any(head["trained"] for head in result["heads"].values())
