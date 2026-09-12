"""Manufactured analytic potentials test software, never research energetics."""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.mep import NEB

from pincer_catmech.kinetics.neb_ts_search import (
    FREQUENCY_FACTOR, full_hessian, refine_dimer, run_path, transfer_mode_overlap,
)


class PairSpringCalculator(Calculator):
    """Exactly differentiated pair springs; negative k supplies a test saddle."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, reference, springs=None):
        super().__init__()
        self.pairs = [(i, j, reference.get_distance(i, j), 1.0)
                      for i in range(len(reference)) for j in range(i)] if springs is None else springs
        self.calls = 0

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        energy, forces = 0.0, np.zeros((len(atoms), 3))
        for i, j, equilibrium, stiffness in self.pairs:
            vector = atoms.positions[i] - atoms.positions[j]
            distance = np.linalg.norm(vector)
            shift = distance - equilibrium
            energy += 0.5 * stiffness * shift**2
            force = -stiffness * shift * vector / distance
            forces[i] += force
            forces[j] -= force
        self.results = {"energy": energy, "forces": forces}


def test_full_hessian_diatomic_analytic_frequency_and_force_count():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.0]])
    atoms.calc = PairSpringCalculator(atoms, [(0, 1, 1.0, 2.0)])
    original = atoms.positions.copy()
    result = full_hessian(atoms, displacement=1e-4)
    reduced_mass = atoms.get_masses()[0] / 2
    assert result.external_rank == 5
    assert len(result.frequencies_cm1) == 1
    assert result.frequencies_cm1[0] == pytest.approx(FREQUENCY_FACTOR * np.sqrt(2 / reduced_mass), rel=1e-8)
    assert atoms.calc.calls == 6 * len(atoms)
    assert np.array_equal(atoms.positions, original)


def test_projection_keeps_every_internal_mode_and_is_rotation_invariant():
    atoms = Atoms("HCO", positions=[[0, 0, 0], [1, 0, 0], [0.4, 1.3, 0]])
    atoms.calc = PairSpringCalculator(atoms)
    result = full_hessian(atoms, 1e-4)
    assert result.external_rank == 6
    assert len(result.frequencies_cm1) == 3 * len(atoms) - 6
    assert np.all(result.frequencies_cm1 > 0)
    moved = atoms.copy()
    moved.rotate(37, [1, 2, 3])
    moved.translate([11, -2, 5])
    moved.calc = PairSpringCalculator(moved)
    rotated = full_hessian(moved, 1e-4)
    assert rotated.frequencies_cm1 == pytest.approx(result.frequencies_cm1, rel=2e-7)


def test_hessian_retains_imaginary_modes_without_absolute_value_conversion():
    atoms = Atoms("HCO", positions=[[0, 0, 0], [1, 0, 0], [.4, 1.3, 0]])
    springs = [(0, 1, atoms.get_distance(0, 1), -2),
               (0, 2, atoms.get_distance(0, 2), 2), (1, 2, atoms.get_distance(1, 2), 2)]
    atoms.calc = PairSpringCalculator(atoms, springs)
    result = full_hessian(atoms, 1e-4)
    assert np.sum(result.frequencies_cm1 < 0) == 1
    assert result.antisymmetry_relative < 1e-7


def test_full_hessian_restores_positions_after_backend_failure():
    class Broken(PairSpringCalculator):
        def calculate(self, *args, **kwargs):
            raise RuntimeError("controlled analytic force failure")

    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    original = atoms.positions.copy()
    atoms.calc = Broken(atoms)
    with pytest.raises(RuntimeError, match="controlled"):
        full_hessian(atoms)
    assert np.array_equal(atoms.positions, original)


def transfer_fixture():
    # Two parallel, deliberately manufactured H-transfer coordinates.
    atoms = Atoms("OHNCHRu", positions=[[-1, 0, 0], [0, 0, 0], [1, 0, 0],
                                         [-1, 3, 0], [0, 3, 0], [1, 3, 0]])
    transfer = dict(proton=1, proton_donor=0, proton_acceptor=2,
                    hydride=4, hydride_donor=3, hydride_acceptor=5)
    return atoms, transfer


def test_dual_transfer_direction_and_wrong_mode_rejection():
    atoms, transfer = transfer_fixture()
    mode = np.zeros((len(atoms), 3))
    mode[[1, 4], 0] = 1
    result = transfer_mode_overlap(atoms, mode, transfer)
    assert result["proton_overlap"] > .5 and result["hydride_overlap"] > .5
    assert result["concerted_direction"]
    mode[4, 0] = -1
    assert not transfer_mode_overlap(atoms, mode, transfer)["concerted_direction"]
    mode[:] = 0
    mode[[1, 4], 2] = 1
    result = transfer_mode_overlap(atoms, mode, transfer)
    assert result["combined_overlap"] == 0


def test_ase_neb_has_seven_internal_images_and_independent_calculators():
    first = Atoms("H2", positions=[[0, 0, 0], [0, 0, .9]])
    last = first.copy()
    last.positions[1, 2] = 1.1
    images = [first] + [first.copy() for _ in range(7)] + [last]
    for image in images:
        image.calc = PairSpringCalculator(first, [(0, 1, 1.0, 2)])
    neb = NEB(images, climb=True, method="improvedtangent")
    neb.interpolate()
    assert len(images) == 9
    assert len({id(image.calc) for image in images}) == 9
    assert neb.get_forces().shape == (7 * len(first), 3)
    assert np.all(np.isfinite(neb.get_forces()))


def test_dimer_requires_true_stationarity(tmp_path):
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    atoms.calc = PairSpringCalculator(atoms, [(0, 1, 1.0, -2)])
    tangent = np.array([[0, 0, -1], [0, 0, 1]], dtype=float)
    assert refine_dimer(atoms, tangent, tmp_path, fmax=0.01, steps=3)
    assert full_hessian(atoms).frequencies_cm1[0] < 0


def test_nonbracketing_endpoints_fail_without_barrier(tmp_path):
    atoms, transfer = transfer_fixture()
    result = run_path(atoms, atoms.copy(), lambda label: PairSpringCalculator(atoms), tmp_path / "failure",
                      transfer, endpoint_steps=0, neb_steps=0, ts_steps=0, max_corrections=0)
    assert not result.accepted
    assert result.status == "failed"
    assert result.forward_electronic_barrier_ev is None
    assert "do not bracket" in " ".join(result.failure_reasons)
    persisted = json.loads((tmp_path / "failure" / "result.json").read_text())
    assert persisted["accepted_ts_energy_ev"] is None
    assert persisted["calculator_evaluations"] > 0


def test_deadline_failure_is_persisted_and_never_accepted(tmp_path):
    atoms, transfer = transfer_fixture()
    def expired(label):
        raise TimeoutError("controlled backend deadline")
    result = run_path(atoms, atoms.copy(), expired, tmp_path, transfer)
    assert result.status == "deadline_exhausted"
    assert not result.accepted
    assert "controlled backend deadline" in result.failure_reasons[0]


def test_endpoint_mapping_and_fresh_calculators_enforced(tmp_path):
    atoms, transfer = transfer_fixture()
    calculator = PairSpringCalculator(atoms)
    result = run_path(atoms, atoms.copy(), lambda label: calculator, tmp_path, transfer, endpoint_steps=0)
    assert "reused a calculator" in " ".join(result.failure_reasons)
    with pytest.raises(ValueError, match="permutation"):
        run_path(atoms, atoms.copy(), lambda label: PairSpringCalculator(atoms), tmp_path / "bad", transfer, atom_mapping=[0] * len(atoms))
    with pytest.raises(ValueError, match="identical"):
        run_path(atoms, Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]]), lambda label: calculator, tmp_path / "bad", transfer)


@pytest.mark.parametrize("options", [{"n_internal_images": 5}, {"neb_fmax": 0}, {"hessian_step": -1},
                                    {"max_corrections": -1}, {"deadline_seconds": 0}, {"minimum_transfer_overlap": 2},
                                    {"ts_initial_image": 0}, {"ts_initial_image": 8}, {"ts_initial_image": True}])
def test_invalid_protocol_rejected(tmp_path, options):
    atoms, transfer = transfer_fixture()
    with pytest.raises(ValueError):
        run_path(atoms, atoms.copy(), lambda label: PairSpringCalculator(atoms), tmp_path, transfer, **options)


class ManufacturedTransferPotential(Calculator):
    """Invariant pair-distance double well for a full software acceptance test.

    This potential is manufactured analytically, is not GFN2-xTB, and has no
    chemical interpretation. Its exact gradient checks the entire solver chain.
    """

    implemented_properties = ["energy", "forces"]
    center = np.array([[-1, 0, 0], [0, 0, 0], [1, 0, 0],
                       [-1, 3, 0], [0, 3, .25], [1, 3, .5]], dtype=float)
    tangent = np.zeros((6, 3))
    tangent[1] = [.5, 0, 0]
    tangent[4] = .5 * (center[5] - center[3]) / np.linalg.norm(center[5] - center[3])

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        q, gradient_q = 0.0, np.zeros((6, 3))
        for neighbor, sign in ((0, 1), (2, -1)):
            vector = atoms.positions[1] - atoms.positions[neighbor]
            distance = np.linalg.norm(vector)
            q += sign * distance
            gradient_q[1] += sign * vector / distance
            gradient_q[neighbor] -= sign * vector / distance
        reference = self.center + q * self.tangent
        energy, derivative_q = (q**2 - .36)**2, 4 * q * (q**2 - .36)
        gradient = np.zeros((6, 3))
        for i in range(6):
            for j in range(i):
                vector = atoms.positions[i] - atoms.positions[j]
                distance = np.linalg.norm(vector)
                ref_vector = reference[i] - reference[j]
                ref_distance = np.linalg.norm(ref_vector)
                residual = distance - ref_distance
                energy += 2.5 * residual**2
                pair_gradient = 5 * residual * vector / distance
                gradient[i] += pair_gradient
                gradient[j] -= pair_gradient
                derivative_q -= 5 * residual * np.dot(ref_vector, self.tangent[i] - self.tangent[j]) / ref_distance
        self.results = {"energy": energy, "forces": -gradient - derivative_q * gradient_q}


def test_complete_neb_dimer_hessian_and_descent_on_manufactured_potential(tmp_path):
    _, transfer = transfer_fixture()
    reactant = Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center - .6 * ManufacturedTransferPotential.tangent)
    product = Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + .6 * ManufacturedTransferPotential.tangent)
    result = run_path(reactant, product, lambda label: ManufacturedTransferPotential(), tmp_path, transfer,
                      endpoint_steps=0, neb_steps=120, ts_steps=80, descent_steps=180,
                      max_corrections=0, hessian_step=.001)
    assert result.accepted, result.failure_reasons
    attempt = result.attempts[-1]
    assert attempt["imaginary_count"] == 1
    assert -1800 <= min(attempt["frequencies_cm1"]) <= -300
    assert attempt["transfer_mode"]["concerted_direction"]
    assert {branch["nearest_endpoint"] for branch in attempt["descent_branches"]} == {0, 1}
    assert result.forward_electronic_barrier_ev == pytest.approx(.1296, abs=.003)
    assert (tmp_path / "attempt_00" / "full_hessian.npz").is_file()


def test_restart_uses_saved_frames_fresh_forces_and_full_acceptance(tmp_path, monkeypatch):
    from ase.io import read, write
    import hashlib
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    source = tmp_path / "previous.traj"
    write(source, frames)
    def forbidden_interpolation(*args, **kwargs):
        raise AssertionError("Restart must preserve supplied coordinates")
    monkeypatch.setattr(NEB, "interpolate", forbidden_interpolation)
    result = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "continued", transfer,
                      restart_images=frames, restart_source=source, endpoint_steps=0, neb_steps=120,
                      ts_steps=80, descent_steps=180, max_corrections=0, hessian_step=.001)
    assert result.accepted, result.failure_reasons
    assert result.endpoint_converged == (True, True) and result.calculator_evaluations > 6 * len(frames[0])
    protocol = json.loads((tmp_path / "continued" / "protocol.json").read_text())
    assert protocol["restart"]["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    restored = read(tmp_path / "continued" / "restart_input.traj", ":")
    assert all(np.array_equal(a.positions, b.positions) for a, b in zip(frames, restored))


def test_restart_rejects_changed_mapping_state_and_endpoints(tmp_path):
    atoms, transfer = transfer_fixture()
    frames = [atoms.copy() for _ in range(9)]
    def invoke(saved):
        return run_path(atoms, atoms.copy(), lambda label: ManufacturedTransferPotential(), tmp_path,
                        transfer, restart_images=saved)
    with pytest.raises(ValueError, match="nine"):
        invoke(frames[:-1])
    frames[4].info["charge"] = 1
    with pytest.raises(ValueError, match="electronic-state"):
        invoke(frames)
    frames[4].info.clear()
    frames[4].numbers[0] = 6
    with pytest.raises(ValueError, match="identities"):
        invoke(frames)
    frames[4] = atoms.copy()
    frames[-1].positions[0, 0] += .01
    with pytest.raises(ValueError, match="endpoint coordinates"):
        invoke(frames)


def test_explicit_unconverged_band_guess_keeps_all_independent_saddle_gates(tmp_path, monkeypatch):
    from ase.optimize import FIRE
    original_run = FIRE.run
    def unfinished_band(self, *args, **kwargs):
        # Only the band convergence flag is forced false. All forces, dimer,
        # full Hessian and downhill optimizations use the analytic calculator.
        if isinstance(self.atoms, NEB):
            return False
        return original_run(self, *args, **kwargs)
    monkeypatch.setattr(FIRE, "run", unfinished_band)
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    kwargs = dict(restart_images=frames, endpoint_steps=0, neb_steps=0, ts_steps=80,
                  descent_steps=180, max_corrections=0, hessian_step=.001)
    rejected = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "default", transfer, **kwargs)
    assert not rejected.accepted and not rejected.neb_converged
    accepted = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "explicit", transfer,
                        refine_unconverged_band=True, **kwargs)
    assert accepted.accepted, accepted.failure_reasons
    assert not accepted.neb_converged and accepted.provisional_neb_initial_guess
    attempt = accepted.attempts[-1]
    assert attempt["initial_guess_from_unconverged_band"] and attempt["imaginary_count"] == 1
    assert attempt["true_fmax_ev_a"] <= .03 and attempt["transfer_mode"]["concerted_direction"]
    assert {branch["nearest_endpoint"] for branch in attempt["descent_branches"]} == {0, 1}


def test_explicit_local_image_remains_distinct_from_global_band_peak(tmp_path):
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    result = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path, transfer,
                      restart_images=frames, ts_initial_image=3, endpoint_steps=0, neb_steps=120,
                      ts_steps=100, descent_steps=180, max_corrections=0, hessian_step=.001)
    assert result.accepted, result.failure_reasons
    assert result.global_peak_image == 4 and result.selected_initial_image == 3
    assert result.image_energies_ev[3] < result.image_energies_ev[4]
    assert result.attempts[-1]["imaginary_count"] == 1


def test_saved_dimer_candidate_preserves_geometry_and_all_saddle_gates(tmp_path):
    from ase.io import read, write
    import hashlib
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    candidate = Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center - .12 * ManufacturedTransferPotential.tangent)
    source = tmp_path / "actual_previous_dimer.traj"
    write(source, candidate)
    output = tmp_path / "continued"
    result = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), output, transfer,
                      restart_images=frames, restart_ts_candidate=candidate, restart_ts_source=source,
                      endpoint_steps=0, neb_steps=120, ts_steps=100, descent_steps=180, max_corrections=0, hessian_step=.001)
    assert result.accepted, result.failure_reasons
    assert result.restarted_ts_candidate and result.selected_initial_image is None
    assert np.array_equal(read(output / "restart_ts_input.traj").positions, candidate.positions)
    protocol = json.loads((output / "protocol.json").read_text())
    assert protocol["restart_ts_candidate"]["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="requires restart_images"):
        run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "invalid", transfer,
                 restart_ts_candidate=candidate, restart_ts_source=source)


@pytest.mark.parametrize("omit_mode", [False, True])
def test_complete_hessian_provider_is_reprojected_and_all_saddle_gates_remain(tmp_path, omit_mode):
    """Injected analytic full Hessian exercises the native-provider software API."""
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    calls = []
    def provider(candidate, attempt_directory):
        calls.append(attempt_directory)
        result = full_hessian(candidate, .001)
        if omit_mode:
            result.frequencies_cm1 = result.frequencies_cm1[1:]
        return result
    result = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path, transfer,
        restart_images=frames, endpoint_steps=0, neb_steps=120, ts_steps=80, descent_steps=180,
        max_corrections=0, hessian_provider=provider)
    assert len(calls) == 1
    if omit_mode:
        assert not result.accepted and "spectrum disagrees" in " ".join(result.failure_reasons)
    else:
        assert result.accepted, result.failure_reasons
        attempt = result.attempts[0]
        assert attempt["imaginary_count"] == 1 and attempt["hessian_displacement_angstrom"] == .001
        assert attempt["transfer_mode"]["concerted_direction"]
        assert {branch["nearest_endpoint"] for branch in attempt["descent_branches"]} == {0, 1}
        assert len(attempt["hessian_sha256"]) == 64


def test_saved_actual_hessian_mode_initializes_dimer_and_retains_all_gates(tmp_path, monkeypatch):
    import hashlib
    import pincer_catmech.kinetics.neb_ts_search as engine
    from ase.io import write
    _, transfer = transfer_fixture()
    frames = [Atoms("OHNCHRu", positions=ManufacturedTransferPotential.center + q * ManufacturedTransferPotential.tangent)
              for q in np.linspace(-.6, .6, 9)]
    candidate = frames[4].copy()
    candidate.calc = ManufacturedTransferPotential()
    source = tmp_path / "actual_software_hessian.npz"
    spectrum = full_hessian(candidate, .001)
    spectrum.save(source, candidate)
    mode = spectrum.cartesian_modes[np.argmin(spectrum.frequencies_cm1)] * 3.7
    candidate_source = tmp_path / "candidate.traj"
    write(candidate_source, candidate)
    initial_modes = []
    original_refine = engine.refine_dimer
    def inspect_initial_mode(atoms, supplied, *args):
        initial_modes.append(supplied.copy())
        return original_refine(atoms, supplied, *args)
    monkeypatch.setattr(engine, "refine_dimer", inspect_initial_mode)
    result = run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "run", transfer,
        restart_images=frames, restart_ts_candidate=candidate, restart_ts_source=candidate_source,
        restart_ts_mode=mode, restart_ts_mode_source=source,
        endpoint_steps=0, neb_steps=120, ts_steps=80, descent_steps=180, max_corrections=0, hessian_step=.001)
    assert result.accepted, result.failure_reasons
    assert np.allclose(initial_modes[0], mode/np.linalg.norm(mode))
    protocol = json.loads((tmp_path / "run/protocol.json").read_text())
    assert protocol["restart_ts_mode"]["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    with np.load(tmp_path / "run/restart_ts_mode.npz") as saved:
        assert np.array_equal(saved["input_mode"], mode)
        assert np.allclose(saved["normalized_mode"], initial_modes[0])
    assert result.attempts[0]["imaginary_count"] == 1
    assert result.attempts[0]["transfer_mode"]["proton_overlap"] >= .15
    assert {branch["nearest_endpoint"] for branch in result.attempts[0]["descent_branches"]} == {0, 1}
    for invalid in (np.zeros((6,3)), np.ones((6,2)), np.full((6,3), np.nan), np.full((6,3), np.inf)):
        with pytest.raises(ValueError, match="restart_ts_mode"):
            run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "invalid", transfer,
                restart_images=frames, restart_ts_candidate=candidate, restart_ts_source=candidate_source,
                restart_ts_mode=invalid, restart_ts_mode_source=source)
    with pytest.raises(ValueError, match="requires restart_ts_candidate"):
        run_path(frames[0], frames[-1], lambda label: ManufacturedTransferPotential(), tmp_path / "missing", transfer,
                 restart_ts_mode=mode, restart_ts_mode_source=source)
