"""Cross-module numerical contract on a physically generated solvent cluster."""
from pathlib import Path
import numpy as np
import torch
from ase.io import read

from pincer_catmech.generators.solvation_clusters import build_solvation_cluster, hydrogen_bond_metrics
from pincer_catmech.kinetics.master_kinetics import MasterKinetics, software_fixture_rates
from pincer_catmech.models.pincer_egnn import PincerEGNN


def test_cluster_e3_symmetry_and_exact_kinetics_jacobian_contract():
    repo = Path(__file__).resolve().parents[1]
    reference = repo / "data/reference_thermochemistry"
    cluster = build_solvation_cluster(3, seed=1,
        substrate=read(reference / "hemiaminal/stationarity/accepted_geometry.xyz"),
        solvent=read(reference / "tert_butanol/stationarity/accepted_geometry.xyz"))
    metrics = hydrogen_bond_metrics(cluster.atoms, cluster.hydrogen_bonds)
    assert all(1.6 <= row["H_acceptor_distance_A"] <= 2.1 for row in metrics)
    x = torch.tensor(cluster.atoms.positions, dtype=torch.float64)
    h = torch.tensor(cluster.atoms.numbers[:, None]/100., dtype=torch.float64)
    torch.manual_seed(301)
    network = PincerEGNN(1, hidden_dim=8, num_layers=2).double()
    angle = .817
    rotation = torch.tensor([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]], dtype=torch.float64)
    before = network(h, x)
    for q in (rotation, -rotation):
        after = network(h, x @ q.T+torch.tensor([3., -8., 2.]))
        assert float((before.scalars-after.scalars).abs().max().detach()) < 1e-6
        torch.testing.assert_close(before.positions @ q.T+torch.tensor([3., -8., 2.]), after.positions, atol=1e-11, rtol=0)
    kinetics = MasterKinetics(software_fixture_rates(), tbuoh_activity=.1, allow_software_fixture=True)
    c = np.linspace(.001, .3, 18)
    step = 1e-6
    finite = np.column_stack([(kinetics.rhs(0, c+step*np.eye(18)[j])-kinetics.rhs(0, c-step*np.eye(18)[j]))/(2*step) for j in range(18)])
    assert np.max(np.abs(kinetics.jacobian(0, c)-finite)) < 1e-5
