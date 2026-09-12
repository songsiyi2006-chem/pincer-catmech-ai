"""Grouped, evidence-gated gradient-boosting ensemble for computed TS barriers.

No training on ground-state descriptors without accepted TS labels. Ligand
backbone/substituent groups stay intact across metals in five-fold validation.
Ensemble disagreement is not a calibrated confidence interval or chemical
accuracy estimate. Pareto cost is an explicit ordinal preference, not money.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from pincer_catmech.kinetics.microkinetics import (
    BALANCES, ComputedFreeEnergy, EvidenceError, KineticResult, _composition_balance,
)

DEFAULT_FEATURES = ("bite_angle_deg", "percent_buried_volume", "d_electron_count", "topological_index")
METAL_PREFERENCE_ORDINAL = {"Mn": 1.0, "Fe": 1.0, "Co": 2.0, "Ru": 3.0}


@dataclass(frozen=True)
class BarrierObservation:
    catalyst_id: str
    backbone: str
    phosphine_substituent: str
    descriptors: Mapping[str, float]
    transition_state: ComputedFreeEnergy
    reactants: tuple[tuple[ComputedFreeEnergy, int], ...]
    step_name: str = "dehydrogenation"
    reference_convention: str = "preassociated_complex"

    @property
    def ligand_group(self) -> str:
        return self.backbone.strip() + "::" + self.phosphine_substituent.strip()

    def validate(self, feature_names: Sequence[str]) -> float:
        if not all(v.strip() for v in (self.catalyst_id, self.backbone, self.phosphine_substituent)):
            raise EvidenceError("Catalyst and ligand-group identity must be explicit")
        if self.step_name != "dehydrogenation" or not self.reactants:
            raise EvidenceError("A dehydrogenation TS and complete stoichiometric reactant reference are required")
        temperature = self.transition_state.temperature
        self.transition_state.validate(temperature, transition_state=True, dehydrogenation=True)
        total = 0.0
        for energy, coefficient in self.reactants:
            if isinstance(coefficient, bool) or not isinstance(coefficient, int) or coefficient < 1:
                raise EvidenceError("Reactant coefficients must be positive integers")
            energy.validate(temperature)
            if energy.protocol != self.transition_state.protocol:
                raise EvidenceError("TS and reactants must use the same method, solvent and reference-state protocol")
            total += coefficient * energy.gibbs_kcal_mol
        if _composition_balance(self.reactants) != _composition_balance([(self.transition_state, 1)]):
            raise EvidenceError("TS label and reactant reference do not have the same elemental composition and charge")
        molecularity = sum(coefficient for _, coefficient in self.reactants)
        if self.reference_convention not in ("preassociated_complex", "separated_reactants") or (
            self.reference_convention == "preassociated_complex" and molecularity != 1
        ) or (self.reference_convention == "separated_reactants" and molecularity < 2):
            raise EvidenceError("Activation-reference convention must match reactant molecularity")
        if any(name not in self.descriptors or not math.isfinite(float(self.descriptors[name])) for name in feature_names):
            raise EvidenceError("All selected descriptors must be present and finite")
        barrier = self.transition_state.gibbs_kcal_mol - total
        if not math.isfinite(barrier) or barrier < 0:
            raise EvidenceError("Negative/nonfinite activation free energy requires scientific review")
        return barrier


def training_readiness(observations: Sequence[BarrierObservation], *,
                       feature_names: Sequence[str] = DEFAULT_FEATURES,
                       minimum_samples: int = 10, folds: int = 5) -> dict:
    """Reject unaccepted labels; report insufficient data without model fitting."""
    if folds != 5:
        raise ValueError("The campaign requires five grouped folds")
    if not isinstance(minimum_samples, int) or minimum_samples < 10:
        raise ValueError("At least ten accepted labels are required")
    if not feature_names or len(set(feature_names)) != len(feature_names):
        raise ValueError("Select unique nonempty descriptor names")
    reasons, labels, valid_indices = [], [], []
    for index, observation in enumerate(observations):
        try:
            labels.append(observation.validate(feature_names))
            valid_indices.append(index)
        except (EvidenceError, ValueError, TypeError) as error:
            reasons.append(f"Observation {index}: {error}")
    groups = {observations[index].ligand_group for index in valid_indices}
    catalysts = {observations[index].catalyst_id for index in valid_indices}
    if len(valid_indices) < minimum_samples or len(catalysts) < minimum_samples:
        reasons.append(f"Need >= {minimum_samples} accepted labels from distinct catalysts; found {len(valid_indices)} labels/{len(catalysts)} catalysts")
    if len(groups) < folds:
        reasons.append(f"Need >= {folds} independent ligand groups; found {len(groups)}")
    if observations:
        methods = {o.transition_state.protocol for o in observations}
        temperatures = {o.transition_state.temperature for o in observations}
        identities = {o.transition_state.source_sha256 for o in observations}
        references = {(o.reference_convention, sum(coefficient for _, coefficient in o.reactants)) for o in observations}
        if len(identities) != len(observations):
            reasons.append("Duplicate TS source labels are not independent training examples")
        if len(methods) != 1 or len(temperatures) != 1:
            reasons.append("A fit requires one computational/solvent/reference-state protocol and one label temperature")
        if len(references) != 1:
            reasons.append("A fit requires one activation-reference convention and reactant molecularity")
        catalyst_groups: dict[str, set[str]] = {}
        for observation in observations:
            catalyst_groups.setdefault(observation.catalyst_id, set()).add(observation.ligand_group)
        if any(len(value) != 1 for value in catalyst_groups.values()):
            reasons.append("One catalyst is assigned inconsistent ligand identities")
    return {"status": "ready" if not reasons else "insufficient_computed_data",
            "reasons": reasons, "accepted_samples": len(valid_indices),
            "independent_catalysts": len(catalysts), "ligand_groups": len(groups)}


@dataclass
class EnsembleSurrogate:
    feature_names: tuple[str, ...]
    estimators: list[Any]
    training_minimum: np.ndarray
    training_maximum: np.ndarray
    temperature: float
    method: str
    evidence_level: str = "surrogate_of_accepted_computed_TS_labels"
    solvent: str | None = None
    solvation_state: str | None = None
    reference_convention: str = "preassociated_complex"

    def predict(self, descriptors: Sequence[Mapping[str, float]]) -> dict:
        if not descriptors:
            return {"mean_kcal_mol": [], "ensemble_std_kcal_mol": [], "outside_training_box": [],
                    "uncertainty_kind": "ensemble_disagreement_not_calibrated"}
        if any(any(name not in row for name in self.feature_names) for row in descriptors):
            raise ValueError("Prediction descriptors are incomplete")
        matrix = np.array([[row[name] for name in self.feature_names] for row in descriptors], dtype=float)
        if not np.all(np.isfinite(matrix)):
            raise ValueError("Prediction descriptors must be finite")
        predictions = np.stack([estimator.predict(matrix) for estimator in self.estimators])
        outside = np.any((matrix < self.training_minimum) | (matrix > self.training_maximum), axis=1)
        return {"mean_kcal_mol": predictions.mean(axis=0).tolist(),
                "ensemble_std_kcal_mol": predictions.std(axis=0, ddof=1).tolist(),
                "outside_training_box": outside.tolist(),
                "uncertainty_kind": "ensemble_disagreement_not_calibrated"}


@dataclass
class TrainingOutcome:
    status: str
    reasons: tuple[str, ...]
    model: EnsembleSurrogate | None = None
    metrics: dict = field(default_factory=dict)
    fold_records: list[dict] = field(default_factory=list)
    observation_provenance: list[dict] = field(default_factory=list)


def _fit_ensemble(matrix: np.ndarray, labels: np.ndarray, groups: np.ndarray,
                  members: int, seed: int) -> list[Any]:
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.pipeline import Pipeline
    random = np.random.default_rng(seed)
    unique = np.unique(groups)
    estimators = []
    for member in range(members):
        sampled_groups = random.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([np.flatnonzero(groups == group) for group in sampled_groups])
        estimator = Pipeline([("regressor", GradientBoostingRegressor(
            n_estimators=120, learning_rate=0.04, max_depth=2,
            min_samples_leaf=2, loss="squared_error", random_state=seed+member,
        ))])
        estimator.fit(matrix[indices], labels[indices])
        estimators.append(estimator)
    return estimators


def train_surrogate(observations: Sequence[BarrierObservation], *,
                    feature_names: Sequence[str] = DEFAULT_FEATURES,
                    minimum_samples: int = 10, ensemble_members: int = 8,
                    seed: int = 2026) -> TrainingOutcome:
    """Five-fold ligand-group CV followed by a final group-bootstrap ensemble.

    No hyperparameter tuning on held-out folds is performed. OOF metrics assess
    the supplied calculated labels only, not measured reaction barriers/yields.
    """
    if not isinstance(ensemble_members, int) or ensemble_members < 3:
        raise ValueError("Uncertainty ensemble requires at least three members")
    if not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    readiness = training_readiness(observations, feature_names=feature_names, minimum_samples=minimum_samples)
    if readiness["status"] != "ready":
        return TrainingOutcome(readiness["status"], tuple(readiness["reasons"]))
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.model_selection import GroupKFold
    matrix = np.array([[row.descriptors[name] for name in feature_names] for row in observations], dtype=float)
    labels = np.array([row.validate(feature_names) for row in observations])
    groups = np.array([row.ligand_group for row in observations])
    predictions, uncertainty = np.zeros(len(labels)), np.zeros(len(labels))
    records = []
    for fold, (train, test) in enumerate(GroupKFold(n_splits=5).split(matrix, labels, groups)):
        train_groups, test_groups = set(groups[train]), set(groups[test])
        if train_groups & test_groups:
            raise RuntimeError("Group leakage detected")
        estimators = _fit_ensemble(matrix[train], labels[train], groups[train], ensemble_members, seed+1000*fold)
        ensemble_predictions = np.stack([estimator.predict(matrix[test]) for estimator in estimators])
        predictions[test] = ensemble_predictions.mean(axis=0)
        uncertainty[test] = ensemble_predictions.std(axis=0, ddof=1)
        records.append({"fold": fold, "train_catalysts": [observations[i].catalyst_id for i in train],
                        "test_catalysts": [observations[i].catalyst_id for i in test],
                        "train_ligand_groups": sorted(train_groups), "test_ligand_groups": sorted(test_groups),
                        "test_mae_kcal_mol": float(mean_absolute_error(labels[test], predictions[test]))})
    metrics = {"grouped_oof_mae_kcal_mol": float(mean_absolute_error(labels, predictions)),
               "grouped_oof_rmse_kcal_mol": float(math.sqrt(mean_squared_error(labels, predictions))),
               "grouped_oof_r2": float(r2_score(labels, predictions)) if np.ptp(labels) > 0 else None,
               "oof_predictions_kcal_mol": predictions.tolist(), "labels_kcal_mol": labels.tolist(),
               "oof_ensemble_std_kcal_mol": uncertainty.tolist(),
               "uncertainty_kind": "group_bootstrap_ensemble_disagreement_not_calibrated",
               "folds": 5, "independent_ligand_groups": len(set(groups)), "n_observations": len(labels)}
    final_estimators = _fit_ensemble(matrix, labels, groups, ensemble_members, seed+10000)
    model = EnsembleSurrogate(tuple(feature_names), final_estimators, matrix.min(axis=0), matrix.max(axis=0),
                              observations[0].transition_state.temperature, observations[0].transition_state.method,
                              solvent=observations[0].transition_state.solvent,
                              solvation_state=observations[0].transition_state.solvation_state,
                              reference_convention=observations[0].reference_convention)
    provenance = [{"catalyst_id": row.catalyst_id, "ligand_group": row.ligand_group,
                   "ts_source_sha256": row.transition_state.source_sha256,
                   "ts_source_path": row.transition_state.source_path,
                   "reactant_source_sha256": [energy.source_sha256 for energy, _ in row.reactants],
                   "transition_state": asdict(row.transition_state),
                   "reference_convention": row.reference_convention,
                   "reactants": [{"coefficient": coefficient, "free_energy": asdict(energy)}
                                 for energy, coefficient in row.reactants]}
                  for row in observations]
    return TrainingOutcome("trained_on_accepted_computed_labels", (), model, metrics, records, provenance)


def save_surrogate(outcome: TrainingOutcome, directory: str | Path) -> tuple[Path, Path]:
    if outcome.model is None or outcome.status != "trained_on_accepted_computed_labels":
        raise EvidenceError("No trained model exists; insufficient data cannot be published as a model")
    import joblib
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    model_path, manifest_path = directory / "barrier_ensemble.joblib", directory / "training_manifest.json"
    joblib.dump(outcome.model, model_path)
    manifest = {"status": outcome.status, "metrics": outcome.metrics, "fold_records": outcome.fold_records,
                "provenance": outcome.observation_provenance, "feature_names": outcome.model.feature_names,
                "temperature": outcome.model.temperature, "method": outcome.model.method,
                "solvent": outcome.model.solvent, "solvation_state": outcome.model.solvation_state,
                "activation_reference_convention": outcome.model.reference_convention,
                "uncertainty": "ensemble disagreement; no calibrated coverage guarantee",
                "target": "computed dehydrogenation activation free energy in kcal/mol"}
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    return model_path, manifest_path


@dataclass(frozen=True)
class CandidatePerformance:
    catalyst_id: str
    metal: str
    kinetics: KineticResult
    temperature: float
    loading_mol_percent: float
    base_equivalents: float


def pareto_front(candidates: Sequence[CandidatePerformance], *,
                 metal_preference_ordinal: Mapping[str, float] = METAL_PREFERENCE_ORDINAL) -> list[dict]:
    """Maximize computed average TOF, minimize supplied ordinal preference score."""
    if not candidates:
        return []
    conditions = {(row.temperature, row.loading_mol_percent, row.base_equivalents,
                   float(row.kinetics.time_seconds[-1]), row.kinetics.initial_concentrations_mol_l,
                   row.kinetics.free_base_concentration_mol_l, row.kinetics.co_shuttle_concentration_mol_l,
                   row.kinetics.computational_protocol) for row in candidates}
    if len(conditions) != 1:
        raise ValueError("Pareto candidates must share temperature/loading/base/time, initial concentrations, shuttle activities and protocol")
    if len({row.catalyst_id for row in candidates}) != len(candidates):
        raise ValueError("Pareto catalyst identifiers must be unique")
    values = []
    for row in candidates:
        if row.kinetics.evidence_level != "computed_accepted_thermochemistry":
            raise EvidenceError("Pareto ranking requires TOF from accepted computed thermochemistry")
        if row.kinetics.computational_protocol is None or len(row.kinetics.initial_concentrations_mol_l) != 11:
            raise EvidenceError("Pareto ranking requires preserved protocol and complete initial-state metadata")
        if row.kinetics.temperature is None or not math.isclose(row.temperature, row.kinetics.temperature, rel_tol=0, abs_tol=1e-6):
            raise EvidenceError("Declared Pareto temperature differs from the actual kinetic calculation")
        inventory = BALANCES @ np.array(row.kinetics.initial_concentrations_mol_l)
        if inventory[1] <= 0 or not math.isclose(row.loading_mol_percent, 100*inventory[0]/inventory[1], rel_tol=1e-8, abs_tol=1e-10):
            raise EvidenceError("Declared Pareto loading differs from the conserved initial catalyst/benzyl inventory")
        if row.kinetics.free_base_concentration_mol_l is None:
            raise EvidenceError("Pareto ranking requires the explicit free-base concentration")
        if row.metal not in metal_preference_ordinal:
            raise ValueError(f"No explicit ordinal preference for metal {row.metal}")
        score, tof = float(metal_preference_ordinal[row.metal]), row.kinetics.average_tof_per_second
        if not math.isfinite(score) or score < 0 or not math.isfinite(tof) or tof < 0:
            raise ValueError("Ordinal score and TOF must be finite and nonnegative")
        values.append((tof, score))
    selected = []
    for i, row in enumerate(candidates):
        tof, score = values[i]
        dominated = any(other_tof >= tof and other_score <= score and
                        (other_tof > tof or other_score < score)
                        for j, (other_tof, other_score) in enumerate(values) if i != j)
        if not dominated:
            selected.append({"catalyst_id": row.catalyst_id, "average_tof_per_second": tof,
                             "metal_preference_ordinal": score,
                             "cost_interpretation": "user-defined ordinal preference, not a market price"})
    return sorted(selected, key=lambda row: (row["metal_preference_ordinal"], -row["average_tof_per_second"]))


__all__ = ["BarrierObservation", "EnsembleSurrogate", "TrainingOutcome", "CandidatePerformance",
           "training_readiness", "train_surrogate", "save_surrogate", "pareto_front"]
