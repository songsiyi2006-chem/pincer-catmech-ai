"""Portable, fail-closed receipts for future certified physical kinetic inputs.

Catalog binding checks exact metadata/geometry bytes, composition and charge.
It is not a full chemical-identity proof: the independent identity, stationarity
and reaction-mode certificates remain mandatory. No calculations are fabricated.
Only summary.json publishes readiness; failed runs cannot replace that pointer.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import uuid

import numpy as np

from pincer_catmech.generators.combinatorial_pincer import BACKBONES, METALS, SUBSTITUENTS
from pincer_catmech.kinetics.master_kinetics import (
    MasterKinetics, MasterThermochemistry, SPECIES, STEP_NAMES, METAL, BASE_EQUIVALENTS,
)
from pincer_catmech.kinetics.microkinetics import ComputedFreeEnergy

EVIDENCE = "computed_conditional_on_mechanism_and_speciation"
BINDING_LIMIT = "Catalog bytes, composition and charge binding is not a full chemical-identity proof."


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))


def _model_sha(entry):
    return hashlib.sha256(_json(entry).encode("utf-8")).hexdigest()


def _read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON field: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8-sig"), object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {value}")))


def _write(path, value):
    Path(path).write_text(_json(value) + "\n", encoding="utf-8")


def _resolve(base, value):
    path = Path(value)
    return (path if path.is_absolute() else base / path).resolve()


def _inside(base, value):
    path = _resolve(base, value)
    if not path.is_relative_to(base.resolve()):
        raise ValueError("Receipt artifact escapes its archive")
    return path


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a number")
    value = float(value)
    if not math.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(f"Invalid {name}")
    return value


def _formula(value):
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", value)
    if not tokens or "".join(symbol + count for symbol, count in tokens) != value:
        raise ValueError("Catalog formula is not a plain elemental formula")
    counts = Counter()
    for symbol, count in tokens:
        counts[symbol] += int(count or 1)
    return dict(counts)


def _binding(entry, base, repo):
    binding = entry["catalyst_binding"]
    catalyst = entry["catalyst_id"]
    catalog = repo / "data/structures" / catalyst / "active"
    metadata_path = _resolve(base, binding["metadata_path"])
    if _sha(metadata_path) != binding["metadata_sha256"] or _sha(catalog / "metadata.json") != binding["metadata_sha256"]:
        raise ValueError("Catalog metadata binding/hash mismatch")
    metadata = _read(metadata_path)
    if metadata.get("catalyst_id") != catalyst or metadata.get("state") != "active":
        raise ValueError("Catalog catalyst/state identity mismatch")
    if f"{metadata.get('metal')}_{metadata.get('backbone')}_{metadata.get('substituent')}" != catalyst:
        raise ValueError("Catalog design labels mismatch")
    active = entry["states"]["cat"]
    expected = _formula(metadata["formula"])
    if dict(active["composition"]) != expected or active["charge"] != metadata["charge"]:
        raise ValueError("Active catalyst composition/charge differs from named catalog design")
    if binding["active_state_source_sha256"] != active["source_sha256"]:
        raise ValueError("Active-state source binding mismatch")
    geometry_path = _resolve(base, binding["geometry_path"])
    if _sha(geometry_path) != binding["geometry_sha256"]:
        raise ValueError("Geometry binding/hash mismatch")
    if (catalog / "best_found.xyz").is_file() and _sha(catalog / "best_found.xyz") != binding["geometry_sha256"]:
        raise ValueError("Geometry differs from named catalog best_found.xyz")
    lines = geometry_path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0])
    if count < 1 or len(lines[2:]) != count:
        raise ValueError("Geometry atom count mismatch")
    composition = Counter()
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4 or not all(math.isfinite(float(v)) for v in fields[1:]):
            raise ValueError("Invalid geometry coordinates")
        composition[fields[0]] += 1
    if dict(composition) != expected:
        raise ValueError("Geometry composition differs from named catalog design")


def _validate(document, base, repo, temperatures, base_grid):
    if not isinstance(document, dict) or not isinstance(document.get("models"), list) or not document["models"]:
        raise ValueError("A nonempty models list is required")
    known = {f"{m}_{b}_{s}" for m in METALS for b in BACKBONES for s in SUBSTITUENTS}
    seen, prepared = set(), []
    for entry in document["models"]:
        catalyst = entry["catalyst_id"]
        if catalyst not in known:
            raise ValueError("Unexpected catalyst label")
        temperature = _number(entry["temperature_K"], "temperature", positive=True)
        total_base = _number(entry["total_added_base_equiv"], "total base")
        if temperature not in temperatures or total_base not in base_grid:
            raise ValueError("Physical model is outside the campaign temperature/base grid")
        key = (catalyst, temperature, total_base)
        if key in seen:
            raise ValueError("Duplicate physical catalyst/temperature/base model")
        seen.add(key)
        if not isinstance(entry["speciation_assumption"], str) or not entry["speciation_assumption"].strip():
            raise ValueError("Explicit free-base/tBuOH speciation assumptions are required")
        activity = _number(entry["tbuoh_activity"], "tBuOH activity")
        initial = np.asarray(entry["initial_concentrations_M"], dtype=float)
        if initial.shape != (18,) or not np.isfinite(initial).all() or np.any(initial < 0) or METAL @ initial <= 0:
            raise ValueError("All 18 finite nonnegative initial concentrations and positive metal are required")
        def records(mapping):
            result = {}
            for name, supplied in mapping.items():
                fields = dict(supplied)
                if fields.get("solvent") != "toluene" or fields.get("solvation_state") != "gsolv":
                    raise ValueError("Campaign certificates must use toluene/gsolv")
                fields["source_path"] = str(_resolve(base, fields["source_path"]))
                fields["composition"] = tuple(tuple(item) for item in fields["composition"])
                fields["imaginary_frequencies_cm1"] = tuple(fields.get("imaginary_frequencies_cm1", ()))
                result[name] = ComputedFreeEnergy(**fields)
            return result
        thermo = MasterThermochemistry(temperature, records(entry["states"]), records(entry["transition_states"]),
                                       entry["irreversible_sink_approximation"])
        rates = thermo.rates()  # validates all 20 states, 18 TSs and barriers before any solver
        _binding(entry, base, repo)
        prepared.append((entry, rates, initial, activity))
    return prepared


def _paths(document):
    for entry in document["models"]:
        for group in ("states", "transition_states"):
            for record in entry[group].values():
                yield record, "source_path", "source_sha256"
        binding = entry["catalyst_binding"]
        yield binding, "metadata_path", "metadata_sha256"
        yield binding, "geometry_path", "geometry_sha256"


def _archive(document, base, run):
    portable = json.loads(_json(document))
    (run / "sources").mkdir()
    manifest = {}
    for record, path_key, hash_key in _paths(portable):
        original = record[path_key]
        raw = _resolve(base, original).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != record[hash_key]:
            raise ValueError("Source changed while archiving")
        relative = f"sources/{digest}.artifact"
        (run / relative).write_bytes(raw)
        manifest[original] = {"path": relative, "sha256": digest, "bytes": len(raw)}
        record[path_key] = relative
    _write(run / "portable_models.json", portable)
    _write(run / "source_manifest.json", manifest)
    return portable


def _metadata(entry, input_sha):
    return {"catalyst_id": entry["catalyst_id"], "temperature_K": entry["temperature_K"],
            "total_added_base_equiv": entry["total_added_base_equiv"], "tbuoh_activity": entry["tbuoh_activity"],
            "speciation_assumption": entry["speciation_assumption"], "initial_concentrations_M": entry["initial_concentrations_M"],
            "input_sha256": input_sha, "model_sha256": _model_sha(entry), "evidence_kind": EVIDENCE,
            "solvent": "toluene", "solvation_state": "gsolv", "points": 500, "time_end_s": 3600.}


def physical_kinetics(input_path, output, repo, temperatures, base_grid):
    """Validate every input before integration, archive exact bytes, atomically publish."""
    source, output, repo = Path(input_path).resolve(), Path(output).resolve(), Path(repo).resolve()
    raw = source.read_bytes()
    input_sha = hashlib.sha256(raw).hexdigest()
    document = _read(source)
    _validate(document, source.parent, repo, temperatures, base_grid)
    physical = output / "kinetics/physical"
    runs = physical / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".pending-", dir=runs))
    run_name = input_sha + "-" + uuid.uuid4().hex[:12]
    final_run = runs / run_name
    (staging / "original_input.json").write_bytes(raw)
    portable = _archive(document, source.parent, staging)
    prepared = _validate(portable, staging, repo, temperatures, base_grid)
    rows = []
    for number, (entry, rates, initial, activity) in enumerate(prepared):
        model = MasterKinetics(rates, tbuoh_activity=activity)
        trajectory = model.integrate(initial, end_time=3600., points=500)
        control = model.control_coefficients(trajectory)
        filename = f"certified_model_{number:04d}.npz"
        file = staging / filename
        metadata = _metadata(entry, input_sha)
        np.savez_compressed(file, time_s=trajectory.time, concentrations_M=trajectory.concentrations,
                            fluxes_M_s=trajectory.fluxes, tof_s=trajectory.tof_per_second,
                            rate_control=control["rate_control"],
                            instantaneous_selectivity_control=control["instantaneous_selectivity_control"],
                            cumulative_product_control=control["cumulative_product_control"],
                            initial_concentrations_M=initial, species=np.array(SPECIES), steps=np.array(STEP_NAMES),
                            metadata_json=np.array(_json(metadata)), metal_error_M=trajectory.metal_error_M,
                            base_equivalent_error_M=trajectory.base_equivalent_error_M)
        rows.append({"catalyst_id": entry["catalyst_id"], "temperature_K": entry["temperature_K"],
                     "total_added_base_equiv": entry["total_added_base_equiv"], "initial_free_base_M": float(initial[17]),
                     "tbuoh_activity": activity, "speciation_assumption": entry["speciation_assumption"],
                     "evidence_kind": EVIDENCE, "terminal_TOF_s-1": float(trajectory.tof_per_second[-1]),
                     "metal_error_M": float(trajectory.metal_error_M), "model_sha256": _model_sha(entry),
                     "initial_concentrations_M": _json(entry["initial_concentrations_M"]),
                     "data_file": (final_run / filename).relative_to(output).as_posix(), "sha256": _sha(file)})
    with (staging / "results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    receipt = {"schema_version": 1, "input_sha256": input_sha, "models": rows,
               "run_directory": final_run.relative_to(output).as_posix(), "binding_limit": BINDING_LIMIT,
               "files": {name: _sha(staging / name) for name in
                         ("original_input.json", "portable_models.json", "source_manifest.json", "results.csv")}}
    _write(staging / "receipt.json", receipt)
    staging.rename(final_run)
    # Independently collect this candidate before changing the published pointer.
    _collect_receipt(receipt, output, repo, temperatures, base_grid)
    temporary = physical / (".summary-" + uuid.uuid4().hex + ".json")
    _write(temporary, receipt)
    os.replace(temporary, physical / "summary.json")
    return rows


def _collect_receipt(receipt, output, repo, temperatures, base_grid):
    if receipt.get("schema_version") != 1:
        raise ValueError("Physical receipt schema is missing or unsupported")
    run = _inside(output, receipt["run_directory"])
    for name in ("original_input.json", "portable_models.json", "source_manifest.json", "results.csv"):
        if _sha(run / name) != receipt["files"][name]:
            raise ValueError(f"Physical archived input/hash mismatch: {name}")
    if _sha(run / "original_input.json") != receipt["input_sha256"]:
        raise ValueError("Original input SHA mismatch")
    original, portable = _read(run / "original_input.json"), _read(run / "portable_models.json")
    manifest = _read(run / "source_manifest.json")
    used = set()
    for record, path_key, hash_key in _paths(original):
        original_path = record[path_key]
        copied = manifest[original_path]
        path = _inside(run, copied["path"])
        if _sha(path) != copied["sha256"] or path.stat().st_size != copied["bytes"] or record[hash_key] != copied["sha256"]:
            raise ValueError("Archived source binding/hash mismatch")
        record[path_key] = copied["path"]
        used.add(original_path)
    if used != set(manifest) or original != portable:
        raise ValueError("Portable model does not preserve original input except source paths")
    prepared = _validate(portable, run, repo, temperatures, base_grid)
    rows = receipt["models"]
    if len(rows) != len(prepared):
        raise ValueError("Physical model receipt count mismatch")
    for row, (entry, rates, initial, activity) in zip(rows, prepared):
        expected = _metadata(entry, receipt["input_sha256"])
        for key in ("catalyst_id", "temperature_K", "total_added_base_equiv", "tbuoh_activity", "speciation_assumption", "evidence_kind", "model_sha256"):
            if row[key] != expected[key]:
                raise ValueError(f"Physical row/input {key} mismatch")
        if row["initial_free_base_M"] != initial[17] or row["initial_concentrations_M"] != _json(entry["initial_concentrations_M"]):
            raise ValueError("Physical row initial concentrations mismatch")
        path = _inside(output, row["data_file"])
        if path.parent != run or _sha(path) != row["sha256"]:
            raise ValueError("Physical NPZ path/hash mismatch")
        with np.load(path, allow_pickle=False) as data:
            if json.loads(str(data["metadata_json"].item())) != expected:
                raise ValueError("Physical NPZ/input metadata mismatch")
            if tuple(data["species"]) != SPECIES or tuple(data["steps"]) != STEP_NAMES:
                raise ValueError("Physical NPZ species/step labels mismatch")
            time, c, flux, tof = (data[k] for k in ("time_s", "concentrations_M", "fluxes_M_s", "tof_s"))
            if time.shape != (500,) or c.shape != (500, 18) or flux.shape != (500, 18) or tof.shape != (500,):
                raise ValueError("Physical NPZ must retain all 500 timepoints")
            if not all(np.isfinite(a).all() for a in (time, c, flux, tof)) or not np.allclose(time, np.linspace(0., 3600., 500), rtol=0, atol=1e-10):
                raise ValueError("Physical trajectory values/time grid mismatch")
            if not np.array_equal(data["initial_concentrations_M"], initial) or not np.allclose(c[0], initial, rtol=0, atol=1e-12):
                raise ValueError("Physical NPZ initial concentrations mismatch")
            for key in ("rate_control", "instantaneous_selectivity_control", "cumulative_product_control"):
                if data[key].shape != (500, 18) or np.isinf(data[key]).any():
                    raise ValueError("Physical sensitivity shape/value mismatch")
            metal_error = float(data["metal_error_M"].item())
            base_error = float(data["base_equivalent_error_M"].item())
            if not 0 <= metal_error < 1e-10 or not 0 <= base_error < 1e-10 or row["metal_error_M"] != metal_error:
                raise ValueError("Physical inventory audit mismatch")
            if c.min() < -1e-10 or np.max(np.abs(c @ METAL - initial @ METAL)) >= 1e-10 or np.max(np.abs(c @ BASE_EQUIVALENTS - initial @ BASE_EQUIVALENTS)) >= 1e-10:
                raise ValueError("Physical archived trajectory violates inventories")
            # Re-evaluate rates on saved concentrations; a rewritten NPZ self-hash
            # cannot relabel a trajectory or invent its TOF independently.
            model = MasterKinetics(rates, tbuoh_activity=activity)
            expected_flux = np.array([model.rates(concentration) for concentration in c])
            if not np.allclose(flux, expected_flux, rtol=1e-10, atol=1e-12) or not np.allclose(tof, flux[:, 5] / (initial @ METAL), rtol=1e-10, atol=1e-12):
                raise ValueError("Physical flux/TOF differs from certified model")
            if row["terminal_TOF_s-1"] != float(tof[-1]):
                raise ValueError("Physical row/NPZ TOF mismatch")
    with (run / "results.csv").open(encoding="utf-8", newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    if csv_rows != [{key: str(value) for key, value in row.items()} for row in rows]:
        raise ValueError("Physical CSV differs from receipt rows")
    return rows


def collect_physical_models(output, repo, temperatures, base_grid):
    """Return only fully revalidated archived inputs, certificates and NPZ rows."""
    output, repo = Path(output).resolve(), Path(repo).resolve()
    receipt_path = output / "kinetics/physical/summary.json"
    if not receipt_path.exists():
        return []
    return _collect_receipt(_read(receipt_path), output, repo, temperatures, base_grid)
