"""Two-surface constrained crossing optimization in Hartree and bohr.

The Harvey et al. (1998), DOI 10.1007/s002140050309, seam-projected
gradient is provided explicitly. Optimization uses a derived equality-
constrained SQP/KKT step, damped BFGS Lagrangian Hessian and an exact L1
plus dynamic quadratic penalty merit function. It is not an implementation
of an unverified 'Chaban-Gordon-Dyall' attribution. A first-order crossing
is not automatically a minimum, and neither is a nonadiabatic rate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import time
from typing import Callable, Protocol

import numpy as np


# A conservative project eligibility screen, not a universal spin-purity test.
SPIN_CONTAMINATION_TOLERANCE = 0.1
_MATCHED_FIELDS = ("symbols", "charge", "method", "basis", "solvation", "solvent",
                   "coordinate_units", "energy_units", "gradient_units")
_PROTOCOL_FIELDS = _MATCHED_FIELDS + (
    "multiplicities", "psi4_version", "scf_e_convergence", "scf_d_convergence",
    "grid_radial", "grid_spherical", "scf_algorithm", "soscf_start_convergence",
    "paired_protocol_sha256",
)


def _finite_scalar(value, name):
    array = np.asarray(value)
    if (array.ndim != 0 or array.dtype.kind not in "fiu"
            or not np.isfinite(array)):
        raise ValueError(f"{name} must be a finite real scalar")
    return float(array)


def _metadata_equal(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_metadata_equal(a[k], b[k]) for k in a)
    return bool(np.array_equal(a, b))


def _validate_quality(metadata):
    for key in ("scf_converged", "quality_pass", "intended_catalyst_identity_pass",
                "diagnostic_label_eligible", "spin_contamination_flag"):
        if key in metadata:
            value = metadata[key]
            expected = key != "spin_contamination_flag"
            if not isinstance(value, (bool, np.bool_)) or bool(value) != expected:
                raise ValueError(f"Surface quality check failed: {key}")
    if metadata.get("status") in {"failed", "error", "unconverged", "not_converged",
                                  "timeout", "budget_exhausted", "scf_failed"}:
        raise ValueError(f"Surface quality check failed: status={metadata['status']}")
    if "gradient_kind" in metadata and metadata["gradient_kind"] != "analytic":
        raise ValueError("Supplied gradients must be analytic")
    for key, expected in (("coordinate_units", "bohr"), ("energy_units", "hartree"),
                          ("gradient_units", "hartree/bohr")):
        if key in metadata and metadata[key] != expected:
            raise ValueError(f"Incompatible surface units: {key}")


def _validate_consistency(reference, candidate):
    if (candidate.state1, candidate.state2) != (reference.state1, reference.state2):
        raise ValueError("Electronic state identities changed during evaluation")
    missing = object()
    for key in _PROTOCOL_FIELDS:
        if not _metadata_equal(reference.metadata.get(key, missing), candidate.metadata.get(key, missing)):
            raise ValueError(f"Surface protocol changed during evaluation: {key}")
    before, after = reference.metadata.get("states"), candidate.metadata.get("states")
    if (before is None) != (after is None):
        raise ValueError("Electronic state provenance disappeared or changed during evaluation")
    if before is not None:
        for a, b in zip(before, after):
            for key in _PROTOCOL_FIELDS + ("multiplicity", "reference", "nalpha", "nbeta", "basis_functions"):
                if not _metadata_equal(a.get(key, missing), b.get(key, missing)):
                    raise ValueError(f"Electronic state protocol changed during evaluation: {key}")


@dataclass
class SurfacePair:
    energy1_hartree: float
    energy2_hartree: float
    gradient1_hartree_bohr: np.ndarray
    gradient2_hartree_bohr: np.ndarray
    state1: str
    state2: str
    spin_dependent: bool
    metadata: dict = field(default_factory=dict)

    def validate(self, shape, *, positions_bohr=None):
        """Reject supplied failed quality/provenance; missing metadata proves nothing.

        Metadata-free analytic surfaces remain usable for mathematical tests.
        The S2 cutoff is a conservative screening policy; passing it does not
        establish wavefunction stability, state tracking, or physical accuracy.
        """
        if self.spin_dependent is not True:
            raise ValueError("A spin-dependent Hamiltonian is required; native GFN2 --uhf is insufficient")
        if (not isinstance(self.state1, str) or not isinstance(self.state2, str)
                or not self.state1.strip() or not self.state2.strip() or self.state1 == self.state2):
            raise ValueError("Two distinct state identities are required")
        for name in ("gradient1_hartree_bohr", "gradient2_hartree_bohr"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ValueError("Each finite analytic gradient must match the input geometry")
            setattr(self, name, value)
        self.energy1_hartree = _finite_scalar(self.energy1_hartree, "energy1_hartree")
        self.energy2_hartree = _finite_scalar(self.energy2_hartree, "energy2_hartree")
        if not isinstance(self.metadata, dict):
            raise ValueError("Surface metadata must be a dictionary")
        _validate_quality(self.metadata)
        states = self.metadata.get("states")
        if states is not None:
            if (not isinstance(states, (list, tuple)) or len(states) != 2
                    or any(not isinstance(state, dict) for state in states)):
                raise ValueError("Surface metadata must describe exactly two states in pair order")
            for i, state in enumerate(states, 1):
                _validate_quality(state)
                multiplicity = state.get("multiplicity")
                expected = state.get("expected_spin_squared")
                if multiplicity is not None:
                    if (isinstance(multiplicity, (bool, np.bool_))
                            or not isinstance(multiplicity, (int, np.integer)) or multiplicity < 1):
                        raise ValueError("State multiplicity must be a positive integer")
                    derived = (int(multiplicity)**2 - 1) / 4
                    if expected is not None and not np.isclose(
                            _finite_scalar(expected, "expected_spin_squared"), derived, atol=1e-12, rtol=0):
                        raise ValueError("Expected S2 disagrees with multiplicity")
                    expected = derived
                    identity = getattr(self, f"state{i}")
                    if identity.startswith("multiplicity_") and identity != f"multiplicity_{multiplicity}":
                        raise ValueError("State identity disagrees with metadata multiplicity")
                if expected is not None:
                    expected = _finite_scalar(expected, "expected_spin_squared")
                    if expected < 0:
                        raise ValueError("Expected S2 must be nonnegative")
                if "spin_squared" in state:
                    s2 = _finite_scalar(state["spin_squared"], "spin_squared")
                    if expected is None or s2 < -1e-8:
                        raise ValueError("Supplied S2 requires a valid pure-state reference")
                    if abs(s2 - expected) > SPIN_CONTAMINATION_TOLERANCE:
                        raise ValueError("Spin contamination exceeds the 0.1 S2 eligibility screen")
                for key, actual in (("energy_hartree", getattr(self, f"energy{i}_hartree")),
                                    ("gradient_hartree_bohr", getattr(self, f"gradient{i}_hartree_bohr"))):
                    if key in state:
                        supplied = np.asarray(state[key], dtype=float)
                        if (supplied.shape != np.shape(actual) or not np.all(np.isfinite(supplied))
                                or not np.allclose(supplied, actual, atol=1e-10, rtol=0)):
                            raise ValueError(f"State metadata disagrees with returned {key}")
            for key in _MATCHED_FIELDS:
                supplied = [m[key] for m in (self.metadata, *states) if key in m]
                if supplied and any(not _metadata_equal(supplied[0], item) for item in supplied[1:]):
                    raise ValueError(f"Electronic surfaces have mismatched {key}")
            if "multiplicities" in self.metadata and not _metadata_equal(
                    self.metadata["multiplicities"], [state.get("multiplicity") for state in states]):
                raise ValueError("Pair and state multiplicities disagree")
        if positions_bohr is not None:
            for metadata in (self.metadata, *(states or ())):
                if "positions_bohr" in metadata:
                    supplied = np.asarray(metadata["positions_bohr"], dtype=float)
                    if (supplied.shape != shape or not np.all(np.isfinite(supplied))
                            or not np.allclose(supplied, positions_bohr, atol=1e-10, rtol=0)):
                        raise ValueError("Surface metadata geometry does not match the evaluated geometry")
        return self

    @property
    def gap(self):
        return float(self.energy1_hartree - self.energy2_hartree)

    @property
    def average_energy(self):
        return float((self.energy1_hartree + self.energy2_hartree) / 2)

    @property
    def average_gradient(self):
        return (self.gradient1_hartree_bohr + self.gradient2_hartree_bohr).ravel() / 2

    @property
    def difference_gradient(self):
        return (self.gradient1_hartree_bohr - self.gradient2_hartree_bohr).ravel()


class SpinSurfaceProvider(Protocol):
    def __call__(self, positions_bohr: np.ndarray) -> SurfacePair: ...


def projected_gradient(pair: SurfacePair, degeneracy_tolerance=1e-10):
    """Return P g_average and lambda*, with P=I-dd^T/(d^T d).

    P g1 = P g2 = P g_average; lambda* = -d.g_average/(d.d).
    No absolute-position dot product enters this translation-invariant test.
    """
    if not np.isfinite(degeneracy_tolerance) or degeneracy_tolerance <= 0:
        raise ValueError("Positive finite degeneracy tolerance required")
    d, g = pair.difference_gradient, pair.average_gradient
    dd = float(d @ d)
    if dd <= degeneracy_tolerance**2:
        raise ValueError("Singular crossing constraint: indistinguishable local gradients")
    multiplier = -float(d @ g) / dd
    return g + multiplier * d, multiplier


def harvey_effective_gradient(pair: SurfacePair, gap_scale=1.0):
    """Projected seam gradient plus scaled gap-restoring component c*d.

    The scale carries inverse-energy units. This effective vector need not
    be the gradient of a scalar potential and is not passed blindly to BFGS.
    """
    if not np.isfinite(gap_scale) or gap_scale <= 0:
        raise ValueError("gap_scale must be positive")
    tangent, _ = projected_gradient(pair)
    return tangent + gap_scale * pair.gap * pair.difference_gradient


def saturated_penalty(pair: SurfacePair, alpha: float, epsilon: float):
    """Complete scalar part of the requested rational penalty, for diagnosis.

    F=A+alpha*c^2/(c^2+epsilon); grad F=gA+
    2*alpha*epsilon*c/(c^2+epsilon)^2*d. alpha is Eh, epsilon Eh^2.
    Its restoring derivative decays for large gaps, so it is not used as the
    primary convergence criterion or attributed to Harvey without evidence.
    """
    if not np.isfinite(alpha) or not np.isfinite(epsilon) or alpha < 0 or epsilon <= 0:
        raise ValueError("Nonnegative finite alpha and positive finite epsilon required")
    gap = pair.gap
    denominator = gap * gap + epsilon
    return (pair.average_energy + alpha * gap * gap / denominator,
            pair.average_gradient + (2 * alpha * epsilon * gap / denominator**2) * pair.difference_gradient)


@dataclass(frozen=True)
class MECPConfig:
    max_iterations: int = 80
    gap_tolerance_hartree: float = 1e-4
    projected_gradient_tolerance: float = 3e-4
    difference_gradient_tolerance: float = 1e-10
    trust_radius_bohr: float = 0.25
    initial_hessian: float = 0.5
    initial_penalty: float = 1.0
    maximum_penalty: float = 1000.0
    max_line_search: int = 12
    max_evaluations: int = 500
    wall_seconds: float = 600.0
    armijo: float = 1e-4

    def __post_init__(self):
        for name in ("max_iterations", "max_line_search", "max_evaluations"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("gap_tolerance_hartree", "projected_gradient_tolerance", "difference_gradient_tolerance",
                     "trust_radius_bohr", "initial_hessian", "initial_penalty", "maximum_penalty", "wall_seconds"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.maximum_penalty < self.initial_penalty or not 0 < self.armijo < 0.5:
            raise ValueError("Invalid penalty range or Armijo parameter")


@dataclass
class MECPResult:
    status: str
    converged: bool
    positions_bohr: np.ndarray
    energy1_hartree: float | None
    energy2_hartree: float | None
    gap_hartree: float | None
    projected_gradient_max: float | None
    lagrange_multiplier: float | None
    difference_gradient_norm: float | None
    iterations: int
    evaluations: int
    wall_seconds: float
    history: list[dict]
    minimum_verified: bool = False
    error: str | None = None
    method: str = "KKT-SQP with damped BFGS and dynamic L1/quadratic merit"
    evidence: str = "first-order crossing search; no spin-orbit or rate calculation"

    def to_dict(self):
        data = asdict(self)
        data["positions_bohr"] = self.positions_bohr.tolist()
        return data


def damped_bfgs(hessian, step, gradient_change):
    """Powell-damped BFGS update preserving positive definiteness."""
    bs = hessian @ step
    sbs = float(step @ bs)
    sy = float(step @ gradient_change)
    if sbs <= 1e-20 or not np.isfinite(sbs + sy):
        return hessian.copy()
    theta = 1.0 if sy >= 0.2 * sbs else 0.8 * sbs / (sbs - sy)
    y = theta * gradient_change + (1 - theta) * bs
    updated = hessian - np.outer(bs, bs) / sbs + np.outer(y, y) / float(step @ y)
    return (updated + updated.T) / 2


def kkt_step(hessian, gradient, difference_gradient, gap):
    """Solve B p + d lambda = -g; d.T p = -gap (one equality)."""
    d = np.asarray(difference_gradient, dtype=float).ravel()
    g = np.asarray(gradient, dtype=float).ravel()
    bg, bd = np.linalg.solve(hessian, np.column_stack((g, d))).T
    denominator = float(d @ bd)
    if denominator <= 1e-20:
        raise ValueError("Rank-deficient KKT crossing constraint")
    multiplier = (float(gap) - float(d @ bg)) / denominator
    return -bg - bd * multiplier, multiplier


def optimize_mecp(provider: SpinSurfaceProvider, positions_bohr,
                  config: MECPConfig | None = None,
                  callback: Callable[[dict], None] | None = None) -> MECPResult:
    """Find a regular first-order constrained crossing; preserve all statuses.

    Backends must keep charge, geometry, Hamiltonian and basis matched, change
    only the electronic state, and reject unconverged SCF gradients. Wall time
    is checked between calls; subprocess backends enforce per-call timeouts.
    """
    cfg = config or MECPConfig()
    x = np.asarray(positions_bohr, dtype=float).copy()
    if x.ndim != 2 or x.shape[1] != 3 or not len(x) or not np.all(np.isfinite(x)):
        raise ValueError("positions_bohr must be a finite nonempty (N,3) array")
    start, evaluations, history, pair = time.monotonic(), 0, [], None
    b = np.eye(x.size) * cfg.initial_hessian
    alpha, rho = cfg.initial_penalty, cfg.initial_penalty
    status, error = "maximum_iterations", None

    def evaluate(at):
        nonlocal evaluations
        if evaluations >= cfg.max_evaluations or time.monotonic() - start > cfg.wall_seconds:
            raise TimeoutError("MECP evaluation or wall-time budget exhausted")
        evaluations += 1
        candidate = provider(at.copy()).validate(at.shape, positions_bohr=at)
        if pair is not None:
            _validate_consistency(pair, candidate)
        return candidate

    def merit(p):
        return p.average_energy + rho * abs(p.gap) + 0.5 * alpha * p.gap**2

    try:
        pair = evaluate(x)
        for iteration in range(cfg.max_iterations + 1):
            d, g, gap = pair.difference_gradient, pair.average_gradient, pair.gap
            if np.linalg.norm(d) <= cfg.difference_gradient_tolerance:
                status = "degenerate_surfaces" if abs(gap) <= cfg.gap_tolerance_hartree else "singular_constraint"
                break
            tangent, multiplier = projected_gradient(pair, cfg.difference_gradient_tolerance)
            record = dict(iteration=iteration, evaluations=evaluations,
                          energy1_hartree=float(pair.energy1_hartree), energy2_hartree=float(pair.energy2_hartree),
                          gap_hartree=gap, projected_gradient_max=float(np.max(np.abs(tangent))),
                          difference_gradient_norm=float(np.linalg.norm(d)), lagrange_multiplier=multiplier,
                          penalty_alpha=alpha, merit_rho=rho)
            history.append(record)
            if callback:
                callback(record.copy())
            if abs(gap) <= cfg.gap_tolerance_hartree and record["projected_gradient_max"] <= cfg.projected_gradient_tolerance:
                status = "converged_first_order_crossing"
                break
            if iteration == cfg.max_iterations:
                break
            p, step_multiplier = kkt_step(b, g, d, gap)
            norm = float(np.linalg.norm(p))
            if norm > cfg.trust_radius_bohr:
                p *= cfg.trust_radius_bohr / norm
            rho = min(cfg.maximum_penalty, max(rho, 1.5 * abs(step_multiplier) + 1))
            # Along an SQP step, d.p reduces c. The one-sided L1 derivative
            # is sign(c)*d.p for c != 0 and abs(d.p) at c == 0.
            def directional():
                dc = float(d @ p)
                return float(g @ p) + rho * (np.sign(gap) * dc if gap else abs(dc)) + alpha * gap * dc
            while directional() >= -1e-16 and alpha < cfg.maximum_penalty:
                alpha = min(cfg.maximum_penalty, alpha * 10)
            slope = directional()
            if slope >= 0:
                status = "merit_not_descent"
                break
            old_merit = merit(pair)
            accepted = None
            step_scale = 1.0
            for _ in range(cfg.max_line_search):
                trial_x = x + step_scale * p.reshape(x.shape)
                trial = evaluate(trial_x)
                if merit(trial) <= old_merit + cfg.armijo * step_scale * slope:
                    accepted = trial_x, trial
                    break
                step_scale *= 0.5
            if accepted is None:
                status = "line_search_failed"
                break
            next_x, next_pair = accepted
            # Same multiplier at both endpoints is essential for a secant of
            # the Lagrangian with respect to geometry, not of lambda itself.
            y = (next_pair.average_gradient - g
                 + step_multiplier * (next_pair.difference_gradient - d))
            b = damped_bfgs(b, (next_x - x).ravel(), y)
            if abs(next_pair.gap) > 0.75 * abs(gap) and abs(next_pair.gap) > cfg.gap_tolerance_hartree:
                alpha = min(cfg.maximum_penalty, alpha * 2)
            x, pair = next_x, next_pair
    except TimeoutError as exc:
        status, error = "budget_exhausted", str(exc)
    except (ValueError, RuntimeError, OSError) as exc:
        status, error = "invalid_surfaces_or_backend_failure", str(exc)
    tangent_max = multiplier = difference_norm = None
    if pair is not None:
        difference_norm = float(np.linalg.norm(pair.difference_gradient))
        if difference_norm > cfg.difference_gradient_tolerance:
            tangent, multiplier = projected_gradient(pair, cfg.difference_gradient_tolerance)
            tangent_max = float(np.max(np.abs(tangent)))
    return MECPResult(status, status == "converged_first_order_crossing", x,
                      None if pair is None else float(pair.energy1_hartree),
                      None if pair is None else float(pair.energy2_hartree),
                      None if pair is None else pair.gap, tangent_max, multiplier,
                      difference_norm, max(0, len(history)-1), evaluations,
                      time.monotonic()-start, history, error=error)


def seam_curvature(provider, positions_bohr, *, displacement_bohr=1e-3,
                   remove_rigid_motion=True, curvature_tolerance=1e-5,
                   gap_tolerance_hartree=1e-4, projected_gradient_tolerance=3e-4,
                   antisymmetry_tolerance=1e-5):
    """Finite-difference constrained Hessian certificate (additional QC calls).

    Project out translations/rotations and d, then evaluate Q.T Hess(L) Q
    by central gradient differences along every column of Q. A certificate
    is separate from convergence; zero/negative modes are reported explicitly.
    Excessive antisymmetry rejects the certificate before symmetrization can
    conceal a nonconservative/noisy gradient field. The absolute antisymmetry
    tolerance has Hartree/bohr^2 units. A single finite-difference step is a
    local numerical check, not a displacement-convergence or noise bound.
    """
    x = np.asarray(positions_bohr, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not len(x) or not np.all(np.isfinite(x)):
        raise ValueError("positions_bohr must be a finite nonempty (N,3) array")
    if displacement_bohr <= 0 or not np.isfinite(displacement_bohr):
        raise ValueError("Positive finite finite-difference displacement required")
    for name, value in (("curvature_tolerance", curvature_tolerance),
                        ("gap_tolerance_hartree", gap_tolerance_hartree),
                        ("projected_gradient_tolerance", projected_gradient_tolerance),
                        ("antisymmetry_tolerance", antisymmetry_tolerance)):
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    pair = provider(x.copy()).validate(x.shape, positions_bohr=x)
    tangent, multiplier = projected_gradient(pair)
    excluded = [pair.difference_gradient]
    if remove_rigid_motion:
        centered = x - np.mean(x, axis=0)
        for axis in np.eye(3):
            excluded.extend([np.tile(axis, (len(x), 1)).ravel(), np.cross(centered, axis).ravel()])
    # Column normalization prevents a large molecular radius or a small gap
    # gradient from numerically erasing the constraint in the rank decision.
    excluded = [vector / np.linalg.norm(vector) for vector in excluded if np.linalg.norm(vector) > 0]
    matrix = np.column_stack(excluded)
    u, singular, _ = np.linalg.svd(matrix, full_matrices=True)
    rank = int(np.sum(singular > 1e-10 * max(1., singular[0])))
    q = u[:, rank:]
    columns = []
    for direction in q.T:
        plus_x, minus_x = x + displacement_bohr * direction.reshape(x.shape), x - displacement_bohr * direction.reshape(x.shape)
        plus = provider(plus_x).validate(x.shape, positions_bohr=plus_x)
        minus = provider(minus_x).validate(x.shape, positions_bohr=minus_x)
        _validate_consistency(pair, plus)
        _validate_consistency(pair, minus)
        lp = plus.average_gradient + multiplier * plus.difference_gradient
        lm = minus.average_gradient + multiplier * minus.difference_gradient
        columns.append(q.T @ (lp - lm) / (2 * displacement_bohr))
    reduced = np.column_stack(columns) if columns else np.empty((0, 0))
    if not np.all(np.isfinite(reduced)):
        raise ValueError("Nonfinite finite-difference seam Hessian")
    antisymmetry = float(np.max(np.abs(reduced-reduced.T))) if reduced.size else 0.
    hessian_symmetry_pass = antisymmetry <= antisymmetry_tolerance
    eigenvalues = np.linalg.eigvalsh((reduced + reduced.T) / 2)
    positive_curvature = bool(eigenvalues.size and np.min(eigenvalues) > curvature_tolerance)
    minimum = bool(positive_curvature and hessian_symmetry_pass and abs(pair.gap) <= gap_tolerance_hartree
                   and np.max(np.abs(tangent)) <= projected_gradient_tolerance)
    return dict(eigenvalues_hartree_bohr2=eigenvalues.tolist(), minimum_verified=minimum,
                positive_curvature=positive_curvature,
                antisymmetry_max=antisymmetry, hessian_symmetry_pass=hessian_symmetry_pass,
                antisymmetry_tolerance_hartree_bohr2=antisymmetry_tolerance,
                tangent_dimension=q.shape[1],
                finite_difference_displacement_bohr=displacement_bohr,
                first_order_gap_hartree=pair.gap, projected_gradient_max=float(np.max(np.abs(tangent))),
                interpretation="single-displacement local curvature check; no noise bound, global minimum or SOC claim")
