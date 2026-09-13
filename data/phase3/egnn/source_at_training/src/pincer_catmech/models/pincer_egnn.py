"""Native PyTorch E(3)-equivariant scalar/coordinate message passing.

Based on Satorras, Hoogeboom & Welling, ICML 2021, equations (3)--(6):
https://proceedings.mlr.press/v139/satorras21a.html . No PyG/DGL dependency.
Coordinates are in angstrom; named supervised scalar targets are in eV.
Forward outputs are raw model values, not evidence of trained chemical heads.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Sequence

import torch
from torch import Tensor, nn

DEFAULT_HEADS = ("activation_barrier_ev", "mecp_gap_ev", "aux_conformer_relative_electronic_energy_ev")


@dataclass
class EGNNOutput:
    scalars: Tensor
    positions: Tensor
    node_features: Tensor
    graph_mask: Tensor
    head_names: tuple[str, ...]

    def head(self, name: str) -> Tensor:
        return self.scalars[:, self.head_names.index(name)]


def _sum_rows(values: Tensor, indices: Tensor, count: int) -> Tensor:
    """Differentiable permutation-invariant segment sum using native index_add."""
    return values.new_zeros((count,) + values.shape[1:]).index_add(0, indices, values)


def complete_edge_index(batch: Tensor, node_mask: Tensor | None = None) -> Tensor:
    """Directed within-graph edges, excluding self loops and masked nodes."""
    if batch.ndim != 1 or batch.dtype != torch.long or torch.any(batch < 0):
        raise ValueError("batch must be a nonnegative one-dimensional int64 tensor")
    mask = torch.ones_like(batch, dtype=torch.bool) if node_mask is None else node_mask
    if mask.shape != batch.shape or mask.dtype != torch.bool:
        raise ValueError("node_mask must be Boolean with one entry per node")
    adjacency = (batch[:, None] == batch[None, :]) & mask[:, None] & mask[None, :]
    adjacency.fill_diagonal_(False)
    return adjacency.nonzero(as_tuple=False).T.contiguous()


class EquivariantMessagePassing(nn.Module):
    """Scalar messages times polar relative vectors; no absolute coordinate MLP.

    The invariant coefficient includes 1/sqrt(1+r²), so every per-layer node
    displacement is bounded by coordinate_step, including coincident nodes.
    Directed edge multiplicity contributes to both aggregation and degree.
    """
    def __init__(self, hidden_dim: int, edge_feature_dim: int = 0,
                 coordinate_step: float = 0.05, distance_scale: float = 2.0):
        super().__init__()
        if hidden_dim < 1 or edge_feature_dim < 0 or not math.isfinite(coordinate_step) or coordinate_step < 0:
            raise ValueError("Invalid message dimensions or coordinate step")
        if not math.isfinite(distance_scale) or distance_scale <= 0:
            raise ValueError("distance_scale must be finite and positive")
        self.coordinate_step, self.distance_scale = coordinate_step, distance_scale
        self.message = nn.Sequential(nn.Linear(2 * hidden_dim + 1 + edge_feature_dim, hidden_dim),
                                     nn.SiLU(), nn.Linear(hidden_dim, hidden_dim), nn.SiLU())
        self.coordinate_weight = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.SiLU(),
                                               nn.Linear(hidden_dim, 1, bias=False))
        nn.init.normal_(self.coordinate_weight[-1].weight, std=1e-3)
        self.node_update = nn.Sequential(nn.Linear(2 * hidden_dim, hidden_dim), nn.SiLU(),
                                         nn.Linear(hidden_dim, hidden_dim))

    def forward(self, hidden: Tensor, positions: Tensor, edge_index: Tensor,
                edge_features: Tensor, node_mask: Tensor) -> tuple[Tensor, Tensor]:
        receiver, sender = edge_index
        relative = positions[receiver] - positions[sender]
        squared_distance = relative.square().sum(dim=-1, keepdim=True)
        if not torch.isfinite(squared_distance).all():
            raise ValueError("Coordinate differences overflowed the selected floating-point dtype")
        radial = torch.log1p(squared_distance / self.distance_scale**2)
        message = self.message(torch.cat((hidden[receiver], hidden[sender], radial, edge_features), dim=-1))
        degree = _sum_rows(torch.ones_like(radial), receiver, len(hidden)).clamp_min(1)
        aggregate = _sum_rows(message, receiver, len(hidden)) / degree
        coefficient = self.coordinate_step * torch.tanh(self.coordinate_weight(message)) / torch.sqrt(1 + squared_distance)
        displacement = _sum_rows(relative * coefficient, receiver, len(hidden)) / degree
        updated_hidden = hidden + self.node_update(torch.cat((hidden, aggregate), dim=-1))
        return updated_hidden * node_mask[:, None], positions + displacement


class PincerEGNN(nn.Module):
    """Packed molecular graphs with independent scalar heads and scalar features.

    Edges are directed (receiver, sender). Supplied edge attributes and graph
    features must themselves be invariant scalars, not raw Cartesian vectors.
    Missing edges mean a complete graph; an explicitly empty edge tensor means
    no interactions. Padding nodes retain their finite input coordinates but
    contribute neither messages nor graph readout. Empty-graph outputs carry
    graph_mask=False and must never become observed zero-valued labels.
    """
    def __init__(self, node_feature_dim: int, edge_feature_dim: int = 0,
                 graph_feature_dim: int = 0, hidden_dim: int = 32, num_layers: int = 3,
                 heads: Sequence[str] = DEFAULT_HEADS, coordinate_step: float = 0.05,
                 distance_scale: float = 2.0):
        super().__init__()
        if node_feature_dim < 1 or edge_feature_dim < 0 or graph_feature_dim < 0 or num_layers < 1:
            raise ValueError("Feature dimensions and layer count are invalid")
        if not heads or len(set(heads)) != len(heads) or any(not name or "." in name for name in heads):
            raise ValueError("Head names must be nonempty, unique ModuleDict-compatible names")
        self.node_feature_dim, self.edge_feature_dim = node_feature_dim, edge_feature_dim
        self.graph_feature_dim, self.head_names = graph_feature_dim, tuple(heads)
        self.embedding = nn.Linear(node_feature_dim, hidden_dim)
        self.layers = nn.ModuleList(EquivariantMessagePassing(hidden_dim, edge_feature_dim,
                                    coordinate_step, distance_scale) for _ in range(num_layers))
        self.readouts = nn.ModuleDict({name: nn.Sequential(nn.Linear(2 * hidden_dim + graph_feature_dim, hidden_dim),
                                      nn.SiLU(), nn.Linear(hidden_dim, 1)) for name in self.head_names})
        self.register_buffer("trained_heads", torch.zeros(len(self.head_names), dtype=torch.bool))
        self.register_buffer("pair_only_heads", torch.zeros(len(self.head_names), dtype=torch.bool))

    def forward(self, node_features: Tensor, positions: Tensor, *, edge_index: Tensor | None = None,
                edge_features: Tensor | None = None, batch: Tensor | None = None,
                node_mask: Tensor | None = None, num_graphs: int | None = None,
                graph_features: Tensor | None = None) -> EGNNOutput:
        n = len(positions)
        if positions.shape != (n, 3) or node_features.shape != (n, self.node_feature_dim):
            raise ValueError("Expected positions[N,3] and node_features[N,F]")
        if not positions.is_floating_point() or positions.dtype != node_features.dtype or positions.device != node_features.device:
            raise ValueError("Positions and node features must share a floating-point dtype and device")
        if not torch.isfinite(positions).all() or not torch.isfinite(node_features).all():
            raise ValueError("Inputs, including finite masked padding, must not contain NaN or infinity")
        batch = torch.zeros(n, dtype=torch.long, device=positions.device) if batch is None else batch
        mask = torch.ones(n, dtype=torch.bool, device=positions.device) if node_mask is None else node_mask
        if batch.shape != (n,) or batch.dtype != torch.long or batch.device != positions.device or torch.any(batch < 0):
            raise ValueError("batch must contain one nonnegative int64 graph index per node")
        if mask.shape != (n,) or mask.dtype != torch.bool or mask.device != positions.device:
            raise ValueError("node_mask must be Boolean with one entry per node")
        inferred = int(batch.max()) + 1 if n else 1
        count = inferred if num_graphs is None else num_graphs
        if isinstance(count, bool) or not isinstance(count, int) or count < 0 or (n and count < inferred):
            raise ValueError("num_graphs cannot omit any indexed graph")
        edges = complete_edge_index(batch, mask) if edge_index is None else edge_index
        if edges.ndim != 2 or edges.shape[0] != 2 or edges.dtype != torch.long or edges.device != positions.device:
            raise ValueError("edge_index must be an int64 [2,E] tensor on the input device")
        if edges.numel() and (torch.any(edges < 0) or torch.any(edges >= n)):
            raise ValueError("An edge references a missing node")
        if edges.numel() and torch.any(batch[edges[0]] != batch[edges[1]]):
            raise ValueError("Cross-graph edges are forbidden")
        attributes = positions.new_zeros((edges.shape[1], self.edge_feature_dim)) if edge_features is None else edge_features
        if attributes.shape != (edges.shape[1], self.edge_feature_dim) or attributes.dtype != positions.dtype or attributes.device != positions.device:
            raise ValueError("edge_features must match the directed edges and scalar feature dimension")
        if not torch.isfinite(attributes).all():
            raise ValueError("Edge scalar features must be finite")
        keep = (edges[0] != edges[1]) & mask[edges[0]] & mask[edges[1]]
        edges, attributes = edges[:, keep], attributes[keep]
        context = positions.new_zeros((count, self.graph_feature_dim)) if graph_features is None else graph_features
        if context.shape != (count, self.graph_feature_dim) or context.dtype != positions.dtype or context.device != positions.device or not torch.isfinite(context).all():
            raise ValueError("graph_features must be finite scalar features with shape [num_graphs,Fg]")
        hidden = self.embedding(node_features) * mask[:, None]
        coordinates = positions
        for layer in self.layers:
            hidden, coordinates = layer(hidden, coordinates, edges, attributes, mask)
        sizes = _sum_rows(mask[:, None].to(hidden.dtype), batch, count)
        summed = _sum_rows(hidden, batch, count)
        pooled = torch.cat((summed, summed / sizes.clamp_min(1), context), dim=-1)
        valid = sizes[:, 0] > 0
        scalars = torch.cat([readout(pooled) for readout in self.readouts.values()], dim=-1)
        return EGNNOutput(scalars * valid[:, None], coordinates, hidden, valid, self.head_names)

    def mark_heads_trained(self, names: Sequence[str], *, pair_only: bool = False) -> None:
        """Called only after supervised fitting; readiness metadata remains required."""
        for name in names:
            self.trained_heads[self.head_names.index(name)] = True
            self.pair_only_heads[self.head_names.index(name)] = pair_only

    @torch.no_grad()
    def predict_trained(self, name: str, *args, **kwargs) -> Tensor:
        if not bool(self.trained_heads[self.head_names.index(name)]):
            raise ValueError(f"Head {name!r} is untrained: no supported supervised labels were fitted")
        if bool(self.pair_only_heads[self.head_names.index(name)]):
            raise ValueError("This head is a pairwise energy potential; use predict_pair_trained with a matched reference")
        output = self(*args, **kwargs)
        return output.head(name).masked_fill(~output.graph_mask, float("nan"))

    @torch.no_grad()
    def predict_pair_trained(self, name: str, candidate_inputs: dict, reference_inputs: dict) -> Tensor:
        """Difference of a shared learned potential, never an absolute energy label.

        The caller must validate matched composition, charge and protocol.
        The campaign dataset builder checks these requirements explicitly.
        """
        index = self.head_names.index(name)
        if not bool(self.trained_heads[index]) or not bool(self.pair_only_heads[index]):
            raise ValueError("A trained pairwise head is required")
        candidate, reference = self(**candidate_inputs), self(**reference_inputs)
        if candidate.scalars.shape != reference.scalars.shape:
            raise ValueError("Candidate and reference batches must have matching graph counts")
        valid = candidate.graph_mask & reference.graph_mask
        return (candidate.head(name) - reference.head(name)).masked_fill(~valid, float("nan"))


def masked_mse_loss(predictions: Tensor, targets: Tensor, observed: Tensor) -> Tensor:
    """Missing targets never enter arithmetic; all-missing batches have zero gradient."""
    if targets.shape != predictions.shape or observed.shape != predictions.shape or observed.dtype != torch.bool:
        raise ValueError("Prediction, target and Boolean observation-mask shapes must agree")
    if not torch.isfinite(targets[observed]).all() or not torch.isfinite(predictions[observed]).all():
        raise ValueError("Observed targets and corresponding predictions must be finite")
    if not observed.any():
        return predictions[observed].sum()
    return (predictions[observed] - targets[observed]).square().mean()


@dataclass
class MaskedTargetNormalizer:
    mean: Tensor
    scale: Tensor
    supported: Tensor
    counts: Tensor

    @classmethod
    def fit(cls, targets: Tensor, observed: Tensor, train_indices: Sequence[int]) -> "MaskedTargetNormalizer":
        if targets.ndim != 2 or observed.shape != targets.shape or observed.dtype != torch.bool or not train_indices:
            raise ValueError("A nonempty explicit training split and a matching observation mask are required")
        indices = torch.as_tensor(train_indices, dtype=torch.long, device=targets.device)
        if indices.unique().numel() != indices.numel() or torch.any(indices < 0) or torch.any(indices >= len(targets)):
            raise ValueError("Training indices must be unique and in range")
        values, mask = targets[indices], observed[indices]
        if not torch.isfinite(values[mask]).all():
            raise ValueError("Observed training labels must be finite")
        count = mask.sum(0)
        mean, scale = targets.new_zeros(targets.shape[1]), targets.new_ones(targets.shape[1])
        for head in range(targets.shape[1]):
            if count[head]:
                column = values[:, head][mask[:, head]]
                mean[head] = column.mean()
                deviation = column.std(unbiased=False)
                scale[head] = deviation if deviation > torch.finfo(targets.dtype).eps else 1
        return cls(mean, scale, count > 0, count)

    def transform(self, targets: Tensor) -> Tensor:
        return (targets - self.mean) / self.scale

    def inverse(self, values: Tensor) -> Tensor:
        restored = values * self.scale + self.mean
        return restored.masked_fill(~self.supported, float("nan"))


def grouped_split(groups: Sequence[str], *, seed: int = 20260913,
                  validation_fraction: float = 1 / 6, test_fraction: float = 1 / 6) -> dict[str, list[int]]:
    """Keep every conformer/state/metal sharing the declared family in one split."""
    if any(not isinstance(group, str) or not group for group in groups):
        raise ValueError("Family identifiers must be nonempty strings")
    unique = sorted(set(groups))
    if len(unique) < 3:
        raise ValueError("At least three nonempty independent families are required")
    if not 0 < validation_fraction < 1 or not 0 < test_fraction < 1 or validation_fraction + test_fraction >= 1:
        raise ValueError("Validation/test fractions must leave an independent training split")
    random.Random(seed).shuffle(unique)
    n_test = max(1, round(len(unique) * test_fraction))
    n_validation = max(1, round(len(unique) * validation_fraction))
    if n_test + n_validation >= len(unique):
        raise ValueError("Requested split leaves no training family")
    assignments = {"test": set(unique[:n_test]), "validation": set(unique[n_test:n_test + n_validation]),
                   "train": set(unique[n_test + n_validation:])}
    return {name: [index for index, group in enumerate(groups) if group in families]
            for name, families in assignments.items()}
