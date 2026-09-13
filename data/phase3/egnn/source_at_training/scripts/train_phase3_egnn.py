"""Evidence-gated native-PyTorch EGNN benchmark on real matched conformers.

This auxiliary task learns E(conformer)-E(reference) for the SAME composition,
charge, nominal spin and ALPB protocol. It never substitutes these labels for
activation free energies or spin gaps. All geometries in a ligand family,
across metals and states, remain together in train/validation/test splitting.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from ase.io import read
import torch

from pincer_catmech.models.pincer_egnn import (
    DEFAULT_HEADS, MaskedTargetNormalizer, PincerEGNN, grouped_split, masked_mse_loss,
)
from pincer_catmech.quantum.identity import catalyst_identity_screen
from pincer_catmech.quantum.xtb_backend import atomic_json, parse_energy

REPO = Path(__file__).resolve().parents[1]
ELEMENTS = (1, 6, 7, 8, 15, 25, 26, 27, 44)
STATES = ("active", "hydrogenated", "protonated_reference")
AUXILIARY = DEFAULT_HEADS[2]


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(content):
    return hashlib.sha256(content).hexdigest()


def json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def write_jsonl(path, records):
    lines = [json_bytes(record) for record in records]
    path.write_bytes(b"\n".join(lines) + (b"\n" if lines else b""))


def verify_native_hashes(geometry, output_bytes, geometry_hash, original_output_hash):
    """Verify original LF stdout hash AND preserve exact on-disk evidence bytes.

    Phase2 hashes subprocess text before Windows Path.write_text translates LF
    to CRLF. Reversing only CRLF is the documented original stdout convention;
    arbitrary whitespace or numbers are never normalized.
    """
    if digest(geometry) != geometry_hash:
        raise ValueError("Geometry differs from its original calculation record SHA")
    if digest(output_bytes.replace(b"\r\n", b"\n")) != original_output_hash:
        raise ValueError("Native stdout differs from its original LF-text calculation SHA")


def node_features(atoms, metadata):
    features = np.zeros((len(atoms), len(ELEMENTS) + 3), dtype=float)
    for index, number in enumerate(atoms.numbers):
        if int(number) not in ELEMENTS:
            raise ValueError(f"Unsupported atomic number in the fixed campaign vocabulary: {number}")
        features[index, ELEMENTS.index(int(number))] = 1
    features[metadata["metal_index"], len(ELEMENTS)] = 1
    features[list(metadata["donor_indices"]), len(ELEMENTS) + 1] = 1
    features[metadata["proton_site_index"], len(ELEMENTS) + 2] = 1
    return features.tolist()


def build_dataset(records_path: Path, output: Path) -> dict:
    """Recheck source bytes and identity, then construct fixed-reference pairs."""
    raw_source = records_path.read_bytes()
    source_hash = digest(raw_source)
    if (output / "dataset_manifest.json").exists():
        existing = json.loads((output / "dataset_manifest.json").read_bytes())
        if existing["records_source_sha256"] != source_hash:
            raise ValueError("Source data changed; use a fresh Phase3 EGNN output directory")
        for filename, expected in existing["dataset_files_sha256"].items():
            if digest((output / filename).read_bytes()) != expected:
                raise ValueError(f"Prepared dataset hash mismatch: {filename}")
        return existing
    output.mkdir(parents=True, exist_ok=True)
    (output / "geometries").mkdir(exist_ok=True)
    (output / "evidence").mkdir(exist_ok=True)
    structures, rejected, seen = [], [], set()
    for line_number, raw_line in enumerate(raw_source.splitlines(), 1):
        if not raw_line.strip():
            continue
        row = json.loads(raw_line)
        identifier = f"{row['catalyst_id']}|{row['state']}|{row['conformer']:04d}"
        reason = None
        try:
            native, metadata = row["native_result"], row["metadata"]
            if not native.get("converged") or not row.get("coordination_retained") or not row.get("chemical_identity", {}).get("retained"):
                raise ValueError("Native convergence, donor coordination or complete identity gate failed")
            if row.get("solvent") != "toluene" or row.get("solvation_state") != "gsolv":
                raise ValueError("Only the common ALPB(toluene)/gsolv protocol is eligible")
            folder = Path(native["output_directory"])
            if not folder.exists():
                folder = REPO.parent / "campaign-runs-alpb-toluene" / row["catalyst_id"] / row["state"] / f"conformer_{row['conformer']:04d}"
            geometry = (folder / "xtbopt.xyz").read_bytes()
            output_bytes = (folder / "xtb.out").read_bytes()
            command_bytes = (folder / "command.json").read_bytes()
            command = json.loads(command_bytes)
            verify_native_hashes(geometry, output_bytes, row["optimized_geometry_sha256"], native["output_sha256"])
            native_text = output_bytes.decode("utf-8", errors="strict")
            if "normal termination of xtb" not in native_text.lower() or command.get("returncode") != 0:
                raise ValueError("Native calculation lacks normal termination")
            argv = command["argv"]
            option = lambda key: argv[argv.index(key) + 1]
            if (option("--gfn") != "2" or float(option("--etemp")) != 300 or float(option("--acc")) != .5
                    or option("--alpb") != "toluene" or argv[argv.index("--alpb") + 2] != "gsolv"):
                raise ValueError("Native electronic Hamiltonian/accuracy/temperature/solvent protocol differs")
            energy = float(native["energy_eV"])
            if not math.isfinite(energy) or abs(parse_energy(native_text) - energy) > 1e-6:
                raise ValueError("Record electronic energy disagrees with the hashed native output")
            atoms = read(folder / "xtbopt.xyz")
            if not np.isfinite(atoms.positions).all() or len(atoms) != native["atoms_count"] or atoms.get_chemical_formula() != metadata["formula"]:
                raise ValueError("Coordinates or composition disagree with original metadata")
            if int(option("--chrg")) != metadata["charge"] or int(option("--uhf")) != metadata["unpaired_electrons"]:
                raise ValueError("Charge/nominal spin differs from its native command")
            identity = catalyst_identity_screen(atoms, metadata)
            if not identity["retained"]:
                raise ValueError("Independent full ligand and ancillary identity recheck failed")
            geometry_key = digest(json_bytes({"Z": atoms.numbers.tolist(), "positions_A": np.round(atoms.positions, 8).tolist()}))
            duplicate_key = (row["catalyst_id"], row["state"], geometry_key)
            if duplicate_key in seen:
                raise ValueError("Duplicate same-state geometry at 1e-8 angstrom coordinate precision")
            seen.add(duplicate_key)
            protocol = {"method": "GFN2-xTB", "solvation": "ALPB(toluene)", "solvation_state": "gsolv",
                        "electronic_temperature_K": 300., "accuracy": .5, "energy_unit": "eV",
                        "spin_scope": "Nominal occupancy only; GFN2 is not a spin-gap label source"}
            protocol_hash = digest(json_bytes(protocol))
            structure_id = digest(json_bytes([identifier, geometry_key]))[:24]
            geometry_file = f"geometries/{structure_id}.xyz"
            (output / geometry_file).write_bytes(geometry)
            native_archive = f"evidence/{structure_id}.xtb.out.gz"
            command_file = f"evidence/{structure_id}.command.json"
            (output / native_archive).write_bytes(gzip.compress(output_bytes, mtime=0))
            (output / command_file).write_bytes(command_bytes)
            structures.append({"id": structure_id, "catalyst_id": row["catalyst_id"], "state": row["state"],
                "conformer": row["conformer"], "family_id": metadata["backbone"] + "|" + metadata["substituent"],
                "composition_group": row["catalyst_id"] + "|" + row["state"] + "|" + atoms.get_chemical_formula()
                                     + f"|q{metadata['charge']}|uhf{metadata['unpaired_electrons']}|{protocol_hash}",
                "charge": metadata["charge"], "unpaired_electrons": metadata["unpaired_electrons"],
                "formula": atoms.get_chemical_formula(), "atomic_numbers": atoms.numbers.tolist(),
                "positions_A": atoms.positions.tolist(), "node_features": node_features(atoms, metadata),
                "graph_features": [metadata["charge"], metadata["unpaired_electrons"] / 4]
                                   + [int(row["state"] == state) for state in STATES],
                "electronic_energy_eV": energy, "protocol": protocol, "protocol_sha256": protocol_hash,
                "geometry": geometry_file, "geometry_sha256": digest(geometry), "coordinate_sha256": geometry_key,
                "identity": identity, "source_record_line": line_number, "source_record_line_sha256": digest(raw_line),
                "native_output_sha256": digest(output_bytes), "native_output_archive": native_archive,
                "original_stdout_lf_sha256": native["output_sha256"],
                "native_archive_sha256": digest((output / native_archive).read_bytes()),
                "native_command": command_file, "native_command_sha256": digest(command_bytes)})
        except (ValueError, KeyError, IndexError, OSError) as error:
            reason = f"{type(error).__name__}: {error}"
        if reason:
            rejected.append({"record": identifier, "source_record_line": line_number,
                             "source_record_line_sha256": digest(raw_line), "reason": reason})
        if line_number % 100 == 0:
            print(json.dumps({"stage": "source_audit", "records_checked": line_number, "eligible_structures": len(structures)}), flush=True)
    groups = defaultdict(list)
    for structure in structures:
        groups[structure["composition_group"]].append(structure)
    pairs, references = [], []
    for composition_group, members in sorted(groups.items()):
        members.sort(key=lambda item: (item["conformer"], item["id"]))
        reference = members[0]  # Chosen by pre-existing index, never by target energy.
        references.append({"composition_group": composition_group, "reference_id": reference["id"],
                           "selection": "Lowest original conformer index, then ID; independent of energy"})
        for candidate in members[1:]:
            assert (candidate["formula"], candidate["charge"], candidate["protocol_sha256"], candidate["atomic_numbers"]) == (
                reference["formula"], reference["charge"], reference["protocol_sha256"], reference["atomic_numbers"])
            target = candidate["electronic_energy_eV"] - reference["electronic_energy_eV"]
            pairs.append({"id": candidate["id"] + "__" + reference["id"], "family_id": candidate["family_id"],
                          "candidate_id": candidate["id"], "reference_id": reference["id"], "composition_group": composition_group,
                          "targets_ev": {DEFAULT_HEADS[0]: None, DEFAULT_HEADS[1]: None, AUXILIARY: target},
                          "label_kind": "auxiliary_matched_conformer_electronic_energy_difference",
                          "label_sources": [{"role": role, "structure_id": item["id"], "geometry_sha256": item["geometry_sha256"],
                                             "record_line_sha256": item["source_record_line_sha256"], "native_output_sha256": item["native_output_sha256"]}
                                            for role, item in (("candidate", candidate), ("reference", reference))]})
    write_jsonl(output / "structures.jsonl", structures)
    write_jsonl(output / "pairs.jsonl", pairs)
    atomic_json(output / "rejected_sources.json", rejected)
    atomic_json(output / "reference_selection.json", references)
    manifest = {"created_utc": utc(), "records_source": str(records_path.resolve()), "records_source_sha256": source_hash,
                "source_record_line_hash_convention": "SHA256 of exact source bytes excluding line-ending bytes",
                "original_stdout_hash_convention": "SHA256 of subprocess LF text; Windows CRLF on-disk bytes separately hashed and archived unchanged",
                "eligible_structures": len(structures), "rejected_records": len(rejected), "auxiliary_pairs": len(pairs),
                "families": dict(Counter(pair["family_id"] for pair in pairs)), "composition_groups": len(groups),
                "reference_policy": "Deterministic lowest original conformer index; self-pairs excluded; no energy-based baseline selection",
                "scientific_scope": "Converged electronic optimizations, not Hessian-certified minima; pairwise auxiliary energies only",
                "unsupported_targets": {DEFAULT_HEADS[0]: "Zero validated TS activation-free-energy labels in Phase2",
                                        DEFAULT_HEADS[1]: "GFN2 energies do not discriminate spin; new independent spin-gap labels required"},
                "dataset_files_sha256": {filename: digest((output / filename).read_bytes()) for filename in
                                         ("structures.jsonl", "pairs.jsonl", "rejected_sources.json", "reference_selection.json")}}
    atomic_json(output / "dataset_manifest.json", manifest)
    return manifest


def pack(structures, ids, *, dtype=torch.float32):
    entries = [structures[identifier] for identifier in ids]
    return {"node_features": torch.tensor([row for entry in entries for row in entry["node_features"]], dtype=dtype),
            "positions": torch.tensor([row for entry in entries for row in entry["positions_A"]], dtype=dtype),
            "batch": torch.tensor([index for index, entry in enumerate(entries) for _ in entry["atomic_numbers"]], dtype=torch.long),
            "num_graphs": len(entries), "graph_features": torch.tensor([entry["graph_features"] for entry in entries], dtype=dtype)}


def pair_predictions(model, structures, pairs, indices):
    chosen = [pairs[index] for index in indices]
    ids = [row["candidate_id"] for row in chosen] + [row["reference_id"] for row in chosen]
    prediction = model(**pack(structures, ids)).scalars
    return prediction[:len(chosen)] - prediction[len(chosen):]


def evaluate(model, structures, pairs, indices, batch_size, training_mean):
    observed, predicted, rows = [], [], []
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            selected = indices[start:start + batch_size]
            values = pair_predictions(model, structures, pairs, selected)[:, 2].tolist()
            for index, value in zip(selected, values):
                target = pairs[index]["targets_ev"][AUXILIARY]
                observed.append(target)
                predicted.append(value)
                rows.append({"pair_id": pairs[index]["id"], "family_id": pairs[index]["family_id"],
                             "target_auxiliary_eV": target, "predicted_auxiliary_eV": value,
                             "activation_barrier_prediction_eV": None, "mecp_gap_prediction_eV": None})
    truth, estimate = np.asarray(observed), np.asarray(predicted)
    metrics = {"samples": len(indices), "mae_eV": float(np.mean(np.abs(truth - estimate))),
               "rmse_eV": float(np.sqrt(np.mean((truth - estimate)**2))),
               "equal_reference_energy_baseline_mae_eV": float(np.mean(np.abs(truth))),
               "training_mean_baseline_mae_eV": float(np.mean(np.abs(truth - training_mean)))}
    return metrics, rows


def run_training(output: Path, *, seed=20260913, epochs=30, batch_size=4, hidden_dim=24,
                 num_layers=2, threads=2, max_seconds=600) -> dict:
    import psutil
    torch.set_num_threads(threads)
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    structures = {row["id"]: row for row in map(json.loads, filter(str.strip, (output / "structures.jsonl").read_text().splitlines()))}
    pairs = list(map(json.loads, filter(str.strip, (output / "pairs.jsonl").read_text().splitlines())))
    if len({pair["family_id"] for pair in pairs}) < 3:
        report = {"status": "insufficient_auxiliary_data", "created_utc": utc(),
                  "eligible_structures": len(structures), "auxiliary_pairs": len(pairs),
                  "heads": {head: {"trained": False} for head in DEFAULT_HEADS},
                  "reason": "At least three distinct catalyst families with observed pairs are required for grouped train/validation/test"}
        atomic_json(output / "readiness.json", report)
        return report
    split = grouped_split([pair["family_id"] for pair in pairs], seed=seed)
    target = torch.tensor([[float("nan") if pair["targets_ev"][head] is None else pair["targets_ev"][head]
                            for head in DEFAULT_HEADS] for pair in pairs], dtype=torch.float32)
    mask = torch.isfinite(target)
    normalizer = MaskedTargetNormalizer.fit(target, mask, split["train"])
    configuration = {"node_feature_dim": len(ELEMENTS) + 3, "graph_feature_dim": 5,
                     "hidden_dim": hidden_dim, "num_layers": num_layers, "heads": list(DEFAULT_HEADS)}
    model = PincerEGNN(**configuration)
    for head in DEFAULT_HEADS[:2]:
        for parameter in model.readouts[head].parameters():
            parameter.requires_grad_(False)
    optimizer = torch.optim.AdamW((parameter for parameter in model.parameters() if parameter.requires_grad), lr=2e-3, weight_decay=1e-5)
    run_dir = output / "runs" / f"seed_{seed}_{int(time.time())}"
    run_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = output / "dataset_manifest.json"
    split_data = {name: {"indices": indices, "pair_ids": [pairs[i]["id"] for i in indices],
                         "family_ids": sorted({pairs[i]["family_id"] for i in indices})} for name, indices in split.items()}
    atomic_json(run_dir / "split.json", split_data)
    process = psutil.Process()
    started, deadline = time.monotonic(), time.monotonic() + max_seconds
    peak_rss = process.memory_info().rss
    history, best_state, best_validation, best_epoch = [], None, float("inf"), None
    baseline_mean = float(normalizer.mean[2])
    normalized_target = normalizer.transform(target)
    partial_epoch = False
    for epoch in range(epochs):
        order = list(split["train"])
        random.Random(seed + epoch).shuffle(order)
        total_loss, total_observed = 0., 0
        model.train()
        for start in range(0, len(order), batch_size):
            if time.monotonic() >= deadline:
                partial_epoch = True
                break
            selected = order[start:start + batch_size]
            prediction = pair_predictions(model, structures, pairs, selected)
            normalized_prediction = normalizer.transform(prediction)
            loss = masked_mse_loss(normalized_prediction, normalized_target[selected], mask[selected])
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10)
            optimizer.step()
            total_loss += float(loss.detach()) * len(selected)
            total_observed += len(selected)
            peak_rss = max(peak_rss, process.memory_info().rss)
            if peak_rss > 500 * 1024**2:
                raise MemoryError("EGNN process exceeded the declared 500 MiB memory budget")
        model.eval()
        validation, _ = evaluate(model, structures, pairs, split["validation"], batch_size, baseline_mean)
        if validation["mae_eV"] < best_validation:
            best_validation, best_epoch = validation["mae_eV"], epoch
            best_state = copy.deepcopy(model.state_dict())
        event = {"epoch": epoch, "training_normalized_mse": total_loss / max(1, total_observed),
                 "observed_training_pairs": total_observed, "validation": validation,
                 "partial_epoch_due_to_time_budget": partial_epoch, "elapsed_seconds": time.monotonic() - started}
        history.append(event)
        atomic_json(run_dir / "learning_history.json", history)
        print(json.dumps(event), flush=True)
        if partial_epoch or epoch - best_epoch >= 8:
            break
    if best_state is None or not any(event["observed_training_pairs"] for event in history):
        raise RuntimeError("No supervised optimizer steps completed")
    model.load_state_dict(best_state)
    model.mark_heads_trained([AUXILIARY], pair_only=True)
    model.eval()
    metrics, predictions = {}, []
    for name, indices in split.items():
        metrics[name], rows = evaluate(model, structures, pairs, indices, batch_size, baseline_mean)
        predictions.extend({**row, "split": name} for row in rows)
    write_jsonl(run_dir / "predictions.jsonl", predictions)
    torch.save({"model_state_dict": model.state_dict(), "configuration": configuration,
                "normalization_mean": normalizer.mean, "normalization_scale": normalizer.scale,
                "normalization_supported": normalizer.supported, "normalization_counts": normalizer.counts,
                "dataset_manifest_sha256": digest(manifest_path.read_bytes()), "seed": seed,
                "prediction_scope": "Only difference of shared auxiliary head between matched candidate/reference graphs is trained"},
               run_dir / "model.pt")
    report = {"status": "auxiliary_benchmark_trained", "created_utc": utc(), "torch_version": torch.__version__,
              "device": "CPU", "threads": threads, "configuration": configuration,
              "parameter_count": sum(parameter.numel() for parameter in model.parameters()), "peak_process_rss_bytes": peak_rss,
              "elapsed_seconds": time.monotonic() - started, "best_validation_epoch": best_epoch,
              "completed_or_partial_epochs": len(history), "split_policy": "Entire backbone|substituent family across metals, states and all conformers stays in one split",
              "normalization": {"fit_split": "train", "mean": normalizer.mean.tolist(), "scale": normalizer.scale.tolist(),
                                "supported": normalizer.supported.tolist(), "counts": normalizer.counts.tolist()},
              "metrics": metrics, "heads": {DEFAULT_HEADS[0]: {"trained": False, "observed_labels": 0, "reason": "No validated TS free-energy barriers"},
                                            DEFAULT_HEADS[1]: {"trained": False, "observed_labels": 0, "reason": "No independent spin-dependent gap training dataset"},
                                            AUXILIARY: {"trained": True, "observed_pairs": len(pairs), "pair_only": True}},
              "scope": "Small-family auxiliary conformer-energy benchmark, not trained activation barriers, spin gaps, catalytic ranking or chemical accuracy validation",
              "limitations": ["Only six backbone/substituent families; one validation and one test family is a weak statistical sample",
                              "Related backbone types can occur in different substituent families; this is not leave-one-backbone-out validation",
                              "Conformer optimizations were not individually Hessian-certified; GFN2 spin occupancy cannot supply spin gaps",
                              "Validation selects the checkpoint; test labels are not used in fitting, normalization or model selection"],
              "dataset_manifest_sha256": digest(manifest_path.read_bytes()),
              "model_sha256": digest((run_dir / "model.pt").read_bytes()),
              "split_sha256": digest((run_dir / "split.json").read_bytes()), "run_directory": str(run_dir.relative_to(output))}
    atomic_json(run_dir / "benchmark.json", report)
    atomic_json(output / "readiness.json", report)
    return report


def summarize_benchmark(output: Path) -> dict:
    """Reload the actual checkpoint and measure geometric identities in float64."""
    report = json.loads((output / "readiness.json").read_bytes())
    manifest = json.loads((output / "dataset_manifest.json").read_bytes())
    summary = {**report, "eligible_structures": manifest["eligible_structures"],
               "rejected_records": manifest["rejected_records"], "auxiliary_pairs": manifest["auxiliary_pairs"],
               "composition_groups": manifest["composition_groups"], "families": manifest["families"]}
    if report["status"] != "auxiliary_benchmark_trained":
        atomic_json(output / "summary.json", summary)
        return summary
    checkpoint_path = output / report["run_directory"] / "model.pt"
    if digest(checkpoint_path.read_bytes()) != report["model_sha256"]:
        raise ValueError("Trained checkpoint SHA changed before symmetry verification")
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    model = PincerEGNN(**checkpoint["configuration"]).double().eval()
    model.load_state_dict(checkpoint["model_state_dict"])
    entries = [json.loads(line) for line in (output / "structures.jsonl").read_text().splitlines() if line.strip()]
    structures = {entry["id"]: entry for entry in entries}
    inputs = pack(structures, [entries[0]["id"], entries[-1]["id"]], dtype=torch.float64)
    torch.manual_seed(20260913)
    angles = (.371, -.893, 1.247)
    a, b, c = (torch.tensor(value, dtype=torch.float64) for value in angles)
    rotation_x = torch.tensor([[1., 0., 0.], [0., a.cos(), -a.sin()], [0., a.sin(), a.cos()]], dtype=torch.float64)
    rotation_y = torch.tensor([[b.cos(), 0., b.sin()], [0., 1., 0.], [-b.sin(), 0., b.cos()]], dtype=torch.float64)
    rotation_z = torch.tensor([[c.cos(), -c.sin(), 0.], [c.sin(), c.cos(), 0.], [0., 0., 1.]], dtype=torch.float64)
    rotation = rotation_z @ rotation_y @ rotation_x
    translation = torch.tensor([3.17, -2.41, .63], dtype=torch.float64)
    errors = {}
    with torch.no_grad():
        original = model(**inputs)
        for name, matrix in (("SO3", rotation), ("O3_improper", -rotation)):
            changed = {**inputs, "positions": inputs["positions"] @ matrix.T + translation}
            transformed = model(**changed)
            errors[name] = {"determinant": float(torch.det(matrix)),
                            "max_scalar_absolute_difference_eV": float((transformed.scalars - original.scalars).abs().max()),
                            "max_coordinate_equivariance_difference_A": float((transformed.positions - (original.positions @ matrix.T + translation)).abs().max())}
        order = torch.randperm(len(inputs["positions"]))
        permuted = model(**{**inputs, "node_features": inputs["node_features"][order],
                            "positions": inputs["positions"][order], "batch": inputs["batch"][order]})
        errors["permutation"] = {"max_scalar_absolute_difference_eV": float((permuted.scalars - original.scalars).abs().max()),
                                  "max_coordinate_equivariance_difference_A": float((permuted.positions - original.positions[order]).abs().max())}
        individual = torch.cat([model(**pack(structures, [identifier], dtype=torch.float64)).scalars
                                for identifier in (entries[0]["id"], entries[-1]["id"])])
        errors["batching"] = {"max_scalar_absolute_difference_eV": float((individual - original.scalars).abs().max())}
    maximum_scalar_error = max(values["max_scalar_absolute_difference_eV"] for values in errors.values())
    maximum_coordinate_error = max(values.get("max_coordinate_equivariance_difference_A", 0.) for values in errors.values())
    verification = {"scope": "Software symmetry of actual trained auxiliary checkpoint; not chemical accuracy or barrier/gap training",
                    "dtype": "float64", "angles_radians": angles, "model_sha256": report["model_sha256"],
                    "structure_ids": [entries[0]["id"], entries[-1]["id"]], "errors": errors,
                    "scalar_tolerance_eV": 1e-6, "max_scalar_error_eV": maximum_scalar_error,
                    "coordinate_tolerance_A": 1e-10, "max_coordinate_error_A": maximum_coordinate_error,
                    "passed": maximum_scalar_error < 1e-6 and maximum_coordinate_error < 1e-10}
    if not verification["passed"]:
        raise AssertionError("Actual EGNN checkpoint failed the declared E(3) symmetry bounds")
    atomic_json(output / "symmetry_verification.json", verification)
    summary["symmetry"] = verification
    test_metrics = report["metrics"]["test"]
    improvement = test_metrics["equal_reference_energy_baseline_mae_eV"] - test_metrics["mae_eV"]
    summary["predictive_assessment"] = {
        "test_mae_improvement_over_equal_reference_eV": improvement,
        "predictive_improvement_demonstrated": improvement > 0,
        "interpretation": "This run did not improve the simple equal-reference-energy baseline; no catalyst screening accuracy is demonstrated"
                          if improvement <= 0 else "A descriptive improvement on one held-out family is not broad chemical validation"}
    summary["finished_utc"] = utc()
    summary["computation_finished"] = True
    summary["split"] = json.loads((output / report["run_directory"] / "split.json").read_bytes())
    summary["split"] = {name: {"pairs": len(values["indices"]), "family_ids": values["family_ids"]}
                        for name, values in summary["split"].items()}
    summary["source_files_sha256"] = {str(path.relative_to(REPO)): digest(path.read_bytes()) for path in
        (Path(__file__), REPO / "src/pincer_catmech/models/pincer_egnn.py", REPO / "tests/test_pincer_egnn.py")}
    atomic_json(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, default=REPO / "data/campaign_alpb_toluene/records_with_identity.jsonl")
    parser.add_argument("--output", type=Path, default=REPO / "data/phase3/egnn")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--hidden-dim", type=int, default=24)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--threads", type=int, choices=(1, 2), default=2)
    parser.add_argument("--max-seconds", type=float, default=600)
    args = parser.parse_args()
    if not 1 <= args.epochs <= 500 or not 1 <= args.batch_size <= 8 or not 1 <= args.hidden_dim <= 64 or not 1 <= args.num_layers <= 4 or not 1 <= args.max_seconds <= 1800:
        parser.error("Use the declared bounded CPU model, batch, epoch and runtime ranges")
    manifest = build_dataset(args.records, args.output)
    print(json.dumps({"dataset_ready": manifest}), flush=True)
    if not args.prepare_only:
        print(json.dumps(run_training(args.output, seed=args.seed, epochs=args.epochs, batch_size=args.batch_size,
                         hidden_dim=args.hidden_dim, num_layers=args.num_layers, threads=args.threads,
                         max_seconds=args.max_seconds)), flush=True)
        print(json.dumps(summarize_benchmark(args.output)), flush=True)


if __name__ == "__main__":
    main()
