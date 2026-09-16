"""Thermodynamic identities, unit boundaries and evidence tamper rejection."""
import hashlib
import math

import pytest

from pincer_catmech.kinetics.solvation_sensitivity import (
    R_KCAL, checked_source, conditional_association, recalculate_g, source_path)


def test_equilibrium_relation_and_stoichiometric_activity_power():
    standard = conditional_association(0.0, 3, 383.15, 1.0)
    diluted = conditional_association(0.0, 3, 383.15, 0.1)
    assert standard["log_activity_ratio_complex_to_bare"] == 0
    assert math.exp(diluted["log_activity_ratio_complex_to_bare"]) == pytest.approx(1e-3)
    assert diluted["effective_association_kcal_mol"] == pytest.approx(3*R_KCAL*383.15*math.log(10))


@pytest.mark.parametrize("dg,n,t,a", [(float("nan"),1,383.15,1), (0,0,383.15,1),
    (0,True,383.15,1), (0,1.5,383.15,1), (0,1,0,1), (0,1,383.15,0), (0,1,383.15,float("inf"))])
def test_invalid_states_rejected(dg,n,t,a):
    with pytest.raises(ValueError):
        conditional_association(dg,n,t,a)


def test_sources_are_relocated_and_tampering_fails(tmp_path):
    folder = tmp_path / "data" / "certified"
    folder.mkdir(parents=True)
    path = folder / "source.npz"
    path.write_bytes(b"a")
    record = {"path": "C:\\old\\checkout\\data\\certified\\source.npz",
              "sha256": hashlib.sha256(b"a").hexdigest(), "bytes": 1}
    assert checked_source(tmp_path, record) == path
    path.write_bytes(b"b")
    with pytest.raises(ValueError, match="hash mismatch"):
        checked_source(tmp_path, record)
    with pytest.raises(ValueError):
        source_path(tmp_path, "data/../outside")


def test_missing_units_fail_closed():
    with pytest.raises(ValueError, match="units"):
        recalculate_g({}, 100)
