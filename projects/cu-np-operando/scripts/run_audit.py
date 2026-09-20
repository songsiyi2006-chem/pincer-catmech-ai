"""Run real software checks and save a machine-readable NON-CHEMICAL audit."""
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import unittest

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from neutral_transport import (InputEvidence, steady_surface, film_flux, plug_flow,
    plug_flow_adaptive, target_ranking_gate, cell_energy, faradaic_efficiency)
from test_neutral_transport import equal_fixture, flow_inputs, exact_equal_outlet


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    results = ROOT / "results"
    results.mkdir(parents=True, exist_ok=True)
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern="test_neutral_transport.py")
    checks = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (results / "neutral_transport_test_log.txt").write_text(stream.getvalue(), encoding="utf-8")
    rates, params = equal_fixture(), flow_inputs()
    exact = exact_equal_outlet(rates, params)
    grids = []
    for cells in [20, 40, 80, 160]:
        run = plug_flow(rates, **params, cells=cells)
        grids.append({"cells": cells, "outlet_a_mol_m3": float(run.a_mol_m3[-1]),
                      "absolute_error_mol_m3": abs(float(run.a_mol_m3[-1]) - exact),
                      "relative_error_to_inlet": abs(float(run.a_mol_m3[-1]) - exact) / params["inlet_a_mol_m3"],
                      "material_balance_residual_mol_m3": run.maximum_material_balance_residual_mol_m3})
    adaptive = plug_flow_adaptive(rates, **params)
    integrated_source_checks = []
    for ca, cp in [(100.0, 5.0), (1.0, 100.0)]:
        local_params = params | {"inlet_a_mol_m3": ca, "inlet_p_mol_m3": cp}
        run = plug_flow(rates, **local_params, cells=80)
        source = [film_flux(rates, a, p, local_params["site_density_mol_m2"],
                  local_params["km_a_m_s"], local_params["km_p_m_s"]).flux_mol_m2_s
                  for a, p in zip(run.a_mol_m3[1:], run.p_mol_m3[1:])]
        integrated = local_params["reactor_volume_m3"] * local_params["area_density_m2_m3"] * sum(source) / 80
        integrated_source_checks.append({"inlet_a_mol_m3": ca, "inlet_p_mol_m3": cp,
             "integrated_surface_source_mol_s": integrated,
             "outlet_net_product_formation_mol_s": run.net_product_formation_mol_s,
             "source_outlet_balance_residual_mol_s": integrated - run.net_product_formation_mol_s})
    surface = steady_surface(rates, 0.1, 0.003)
    film = film_flux(rates, 100, 3, 1e-5, 1e-6, 1e-6)
    electrical = cell_energy([0, 1800, 3600], [3, 3, 3], [2, 2, 2], 0.001)
    evidence = InputEvidence("HYPOTHESIS", "synthetic-neutral-fixture-v1", "neutral_fixture", False, "unit_fixture")
    files = sorted([*ROOT.glob("src/**/*.py"), *ROOT.glob("tests/*.py"), *ROOT.glob("scripts/*.py")])
    audit = {
        "schema_version": 1,
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "claim_kind": "SOFTWARE_VERIFICATION",
        "physical_target_prediction": False,
        "experiment_performed": False,
        "fixture_notice": "All energetic and transport fixture values are arbitrary neutral software inputs, not Cu parameters, measurements, computed chemical energies or target predictions.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "platform": platform.platform()},
        "tests": {"run": checks.testsRun, "failures": len(checks.failures),
                  "errors": len(checks.errors), "skipped": len(checks.skipped),
                  "successful": checks.wasSuccessful()},
        "file_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in files},
        "thermodynamic_cycle_log_residual": rates.log_cycle_closure_residual,
        "surface_site_balance_residual": surface.site_balance_residual,
        "surface_steady_residual_s": surface.steady_residual_s,
        "film_balance_residual_mol_m2_s": film.flux_balance_residual_mol_m2_s,
        "independent_integrated_source_checks": integrated_source_checks,
        "grid_verification": {"analytic_outlet_a_mol_m3": exact,
                              "grids": grids,
                              "observed_error_ratios": [grids[n]["absolute_error_mol_m3"] / grids[n+1]["absolute_error_mol_m3"] for n in range(3)],
                              "adaptive_radau_outlet_a_mol_m3": float(adaptive.a_mol_m3[-1]),
                              "adaptive_absolute_error_mol_m3": abs(float(adaptive.a_mol_m3[-1]) - exact)},
        "separate_electrical_unit_fixture": {**asdict(electrical),
           "formed_product_mol": 0.01, "independently_assumed_electrons_per_product": 2,
           "faradaic_efficiency": faradaic_efficiency(0.01, 7200, 2),
           "linked_to_neutral_cycle": False,
           "additional_boundary": "Excludes pumps, separation, other auxiliaries; full-cell electrical input only."},
        "evidence_gate": target_ranking_gate(evidence),
        "not_implemented": ["Cu reaction stoichiometry", "electrode potential/current coupling",
          "constant-potential free energies", "competitive chemistry", "ionic membrane transport",
          "Nernst-Planck migration", "ohmic conduction", "deactivation", "TEA prices", "physical catalyst ranking"],
    }
    (results / "neutral_transport_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(stream.getvalue())
    print(json.dumps({"tests": audit["tests"], "audit_path": str(results / "neutral_transport_audit.json"),
                      "physical_target_prediction": False}, indent=2))
    return 0 if checks.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
