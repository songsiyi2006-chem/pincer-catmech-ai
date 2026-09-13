"""Receipt tests use algebra-only fake certificates and a stubbed solver.

All artifacts live in pytest temporary directories. Passing these tests is
software verification, never a quantum calculation or a chemical certificate.
"""
from collections import Counter
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pincer_catmech.kinetics.master_kinetics import REACTIONS, STEP_NAMES, MasterKinetics, METAL

SPEC = importlib.util.spec_from_file_location("physical_receipts_test_module", Path(__file__).resolve().parents[1] / "scripts/phase3_physical_evidence.py")
physical = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(physical)
TEMPERATURES, BASES = (400., 410.), (.05, .1)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    repo, inputs, output = (tmp_path / name for name in ("repo", "inputs", "output"))
    inputs.mkdir()
    catalyst = "Ru_macho_pnp_Ph"
    catalog = repo / "data/structures" / catalyst / "active"
    catalog.mkdir(parents=True)
    metadata = {"catalyst_id": catalyst, "metal": "Ru", "backbone": "macho_pnp", "substituent": "Ph",
                "state": "active", "formula": "C10H20NOP2Ru", "charge": 0,
                "test_only": "Algebra fixture; this is not a real catalyst design"}
    physical._write(catalog / "metadata.json", metadata)
    cat = dict(Ru=1, C=10, H=20, N=1, P=2, O=1)
    atoms = [symbol + " 0 0 0" for symbol, count in cat.items() for _ in range(count)]
    (catalog / "best_found.xyz").write_text(str(len(atoms)) + "\nAlgebra-only, no physical geometry\n" + "\n".join(atoms) + "\n", encoding="utf-8")
    def add(*parts):
        result = Counter()
        for part in parts:
            result.update(part)
        return dict(result)
    alcohol, aldehyde = dict(C=7, H=8, O=1), dict(C=7, H=6, O=1)
    amine, imine, product = dict(C=6, H=7, N=1), dict(C=13, H=11, N=1), dict(C=13, H=13, N=1)
    h2 = add(cat, dict(H=2))
    formulas = {"cat": cat, "cat_alcohol": add(cat, alcohol), "cat_h2_aldehyde": add(h2, aldehyde),
                "cat_h2": h2, "cat_h2_imine": add(h2, imine), "cat_product": add(cat, product),
                "dimer": add(h2, h2), "cat_co": add(cat, dict(C=1, O=1)),
                "cat_degraded": dict(Ru=1, C=9, H=17, N=1, P=2, O=1),
                "alcohol": alcohol, "aldehyde": aldehyde, "amine": amine, "hemiaminal": add(aldehyde, amine),
                "imine": imine, "product": product, "water": dict(H=2, O=1), "benzene": dict(C=6, H=6),
                "free_base": dict(C=4, H=9, O=1), "tbuoh": dict(C=4, H=10, O=1), "cleavage_fragment": dict(C=5, H=12, O=1)}
    source = inputs / "algebra_only_certificate.txt"
    source.write_text("FAKE CERTIFICATE FOR RECEIPT UNIT TESTS ONLY. No quantum calculation.", encoding="utf-8")
    def record(composition, charge=0, ts=False):
        return dict(gibbs_kcal_mol=24. if ts else 0., temperature=400., source_path=source.name,
                    source_sha256=physical._sha(source), method="ALGEBRA TEST ONLY", accepted=True,
                    stationarity_verified=True, imaginary_frequencies_cm1=[-500.] if ts else [],
                    connectivity_verified=True, evidence_kind="computed", standard_state="1M", solvent="toluene",
                    solvation_state="gsolv", chemical_identity_verified=True, reaction_mode_verified=True,
                    composition=list(map(list, composition.items())), charge=charge)
    states = {name: record(comp, -1 if name in {"free_base", "cat_degraded"} else 0) for name, comp in formulas.items()}
    transitions = {}
    for name, (left, _) in zip(STEP_NAMES, REACTIONS):
        comp = Counter()
        charge = 0
        for species, count in left.items():
            comp.update({symbol: number * count for symbol, number in formulas[species].items()})
            charge += count * states[species]["charge"]
        transitions[name] = record(dict(comp), charge, True)
    initial = [0.] * 18
    initial[0], initial[9], initial[11], initial[17] = .01, 1., 1., .0123
    entry = dict(catalyst_id=catalyst, temperature_K=400., total_added_base_equiv=.05,
                 tbuoh_activity=.123, speciation_assumption="Algebra unit test, not measured speciation",
                 initial_concentrations_M=initial, states=states, transition_states=transitions,
                 irreversible_sink_approximation="Algebra test only, not physical evidence",
                 catalyst_binding=dict(metadata_path=str(catalog / "metadata.json"), metadata_sha256=physical._sha(catalog / "metadata.json"),
                                       geometry_path=str(catalog / "best_found.xyz"), geometry_sha256=physical._sha(catalog / "best_found.xyz"),
                                       active_state_source_sha256=physical._sha(source)))
    document = {"test_only": "Algebra fixtures, no real physical calculations", "models": [entry]}
    input_path = inputs / "input.json"
    physical._write(input_path, document)
    calls = []
    class StubSolver(MasterKinetics):
        def integrate(self, initial, *, end_time, points):
            calls.append("solver")
            c = np.tile(initial, (points, 1))
            flux = np.array([self.rates(row) for row in c])
            return SimpleNamespace(time=np.linspace(0., end_time, points), concentrations=c,
                                   fluxes=flux, tof_per_second=flux[:, 5] / (np.array(initial) @ METAL),
                                   metal_error_M=0., base_equivalent_error_M=0.)
        def control_coefficients(self, trajectory):
            return {key: np.zeros((500, 18)) for key in ("rate_control", "instantaneous_selectivity_control", "cumulative_product_control")}
    monkeypatch.setattr(physical, "MasterKinetics", StubSolver)
    return SimpleNamespace(repo=repo, inputs=inputs, output=output, source=source, path=input_path,
                           document=document, entry=entry, calls=calls, catalog=catalog)


def run(f):
    physical._write(f.path, f.document)
    return physical.physical_kinetics(f.path, f.output, f.repo, TEMPERATURES, BASES)


def collect(f):
    return physical.collect_physical_models(f.output, f.repo, TEMPERATURES, BASES)


def test_portable_original_bytes_and_500_outputs(fixture):
    f = fixture
    rows = run(f)
    receipt = physical._read(f.output / "kinetics/physical/summary.json")
    archive = f.output / receipt["run_directory"]
    assert (archive / "original_input.json").read_bytes() == f.path.read_bytes()
    assert receipt["binding_limit"] == physical.BINDING_LIMIT
    f.source.rename(f.source.with_suffix(".moved"))
    f.path.rename(f.path.with_suffix(".moved"))
    assert collect(f) == rows
    assert f.calls == ["solver"]
    with np.load(f.output / rows[0]["data_file"], allow_pickle=False) as data:
        assert data["concentrations_M"].shape == (500, 18)
        assert data["initial_concentrations_M"][17] == .0123
        assert json.loads(data["metadata_json"].item())["model_sha256"] == rows[0]["model_sha256"]


@pytest.mark.parametrize("change", [
    lambda e: e.update(catalyst_id="Ru_not_a_design_Ph"),
    lambda e: e.update(temperature_K=401.),
    lambda e: e.update(total_added_base_equiv=.051),
    lambda e: e.update(tbuoh_activity=float("inf")),
    lambda e: e.update(speciation_assumption=""),
    lambda e: e.update(initial_concentrations_M=[0.] * 18),
    lambda e: e["states"]["cat"].update(solvent="water"),
    lambda e: e["states"]["cat"].update(solvation_state="gas"),
    lambda e: e["states"]["cat"].update(temperature=410.),
    lambda e: e["states"]["cat"].update(accepted=False),
    lambda e: e["transition_states"].pop(STEP_NAMES[-1]),
    lambda e: e["transition_states"][STEP_NAMES[0]].update(connectivity_verified=False),
    lambda e: e["catalyst_binding"].update(active_state_source_sha256="0" * 64),
    lambda e: e["catalyst_binding"].update(metadata_sha256="0" * 64),
    lambda e: e["catalyst_binding"].update(geometry_sha256="0" * 64),
])
def test_invalid_inputs_fail_before_any_solver(fixture, change):
    f = fixture
    change(f.entry)
    with pytest.raises((ValueError, KeyError)):
        run(f)
    assert not f.calls
    assert not (f.output / "kinetics/physical/summary.json").exists()


def test_duplicate_grid_and_invalid_later_entry_do_not_publish(fixture):
    f = fixture
    f.document["models"].append(copy.deepcopy(f.entry))
    with pytest.raises(ValueError, match="Duplicate physical"):
        run(f)
    assert not f.calls
    f.document["models"][1]["total_added_base_equiv"] = .1
    f.document["models"][1]["transition_states"] = {}
    with pytest.raises(ValueError, match="20 state"):
        run(f)
    assert not f.calls


def test_prior_ready_receipt_survives_failed_solver(fixture, monkeypatch):
    f = fixture
    original = run(f)
    receipt_path = f.output / "kinetics/physical/summary.json"
    raw = receipt_path.read_bytes()
    def broken(*args, **kwargs):
        raise RuntimeError("intentional algebra solver failure")
    monkeypatch.setattr(physical.MasterKinetics, "integrate", broken)
    with pytest.raises(RuntimeError, match="intentional"):
        run(f)
    assert receipt_path.read_bytes() == raw
    assert collect(f) == original


@pytest.mark.parametrize("key,value", [("catalyst_id", "Mn_macho_pnp_Ph"), ("temperature_K", 410.),
    ("total_added_base_equiv", .1), ("tbuoh_activity", .999), ("speciation_assumption", "changed"),
    ("initial_free_base_M", .05), ("terminal_TOF_s-1", 123.)])
def test_collector_rejects_relabelled_rows(fixture, key, value):
    f = fixture
    run(f)
    path = f.output / "kinetics/physical/summary.json"
    receipt = physical._read(path)
    receipt["models"][0][key] = value
    physical._write(path, receipt)
    with pytest.raises(ValueError, match="mismatch"):
        collect(f)


@pytest.mark.parametrize("field", ["metadata_json", "tof_s", "initial_concentrations_M", "fluxes_M_s"])
def test_collector_rejects_npz_tamper_even_with_updated_selfhash(fixture, field):
    f = fixture
    run(f)
    pointer = f.output / "kinetics/physical/summary.json"
    receipt = physical._read(pointer)
    file = f.output / receipt["models"][0]["data_file"]
    with np.load(file, allow_pickle=False) as data:
        arrays = {name: data[name].copy() for name in data.files}
    if field == "metadata_json":
        meta = json.loads(arrays[field].item())
        meta["tbuoh_activity"] = .999
        arrays[field] = np.array(physical._json(meta))
    else:
        arrays[field].flat[-1] += 1.
    np.savez_compressed(file, **arrays)
    receipt["models"][0]["sha256"] = physical._sha(file)
    physical._write(pointer, receipt)
    with pytest.raises(ValueError):
        collect(f)


def test_collector_revalidates_archived_certificate_bytes(fixture):
    f = fixture
    run(f)
    receipt = physical._read(f.output / "kinetics/physical/summary.json")
    archive = f.output / receipt["run_directory"]
    models = physical._read(archive / "portable_models.json")
    source = archive / models["models"][0]["states"]["cat"]["source_path"]
    source.write_text("tampered", encoding="utf-8")
    with pytest.raises(ValueError, match="source binding/hash"):
        collect(f)


def test_no_receipt_returns_no_physical_predictions(fixture):
    assert collect(fixture) == []
