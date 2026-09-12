"""Analytical software fixtures only; these are never campaign energy data."""
from pathlib import Path
import subprocess

import numpy as np
import pytest
from ase import Atoms
from ase.units import Bohr, Hartree

from pincer_catmech.quantum import xtb_backend as backend


def test_energy_last_printed_value_and_fortran_exponent():
    assert backend.parse_energy("TOTAL ENERGY -1.25 Eh\nTOTAL ENERGY -2.0d+00 Eh") == -2 * Hartree


@pytest.mark.parametrize("text", ["", "TOTAL ENERGY NaN Eh", "TOTAL ENERGY -1.0 eV"])
def test_missing_or_invalid_native_energy_fails(text):
    with pytest.raises(ValueError):
        backend.parse_energy(text)


def test_gradient_force_sign_units_and_final_block(tmp_path):
    path = tmp_path / "gradient"
    path.write_text("$grad\n cycle 1\n0 0 0 H\n1.0D-2 2.0D-2 -3.0D-2\n$end\n")
    assert np.allclose(backend.parse_gradient(path, 1), -np.array([[.01, .02, -.03]]) * Hartree / Bohr)


@pytest.mark.parametrize("charge,unpaired", [(0, 1), (True, 0), (0.0, 0), (0, -1)])
def test_electron_occupations_do_not_silently_change(charge, unpaired):
    with pytest.raises(ValueError):
        backend.validate_electrons(Atoms("H2", positions=[[0, 0, 0], [0, 0, .74]]), charge, unpaired)


def test_previous_native_output_is_never_replaced(tmp_path):
    output = tmp_path / "xtb.out"
    output.write_text("earlier scientific evidence")
    with pytest.raises(FileExistsError):
        backend.optimize_geometry(Atoms("He"), tmp_path, executable="unused.exe")
    assert output.read_text() == "earlier scientific evidence"


def test_timeout_preserves_partial_native_output(tmp_path, monkeypatch):
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1, output=b"partial native stdout", stderr=b"partial stderr")
    monkeypatch.setattr(backend.subprocess, "run", timeout)
    with pytest.raises(subprocess.TimeoutExpired):
        backend._execute(Atoms("He"), tmp_path, "unused.exe", 0, 0, 1, ["--grad"], 1)
    assert "partial native stdout" in (tmp_path / "xtb.out").read_text()
    assert "partial stderr" in (tmp_path / "xtb.out").read_text()


def test_restart_is_explicit_and_hash_recorded(tmp_path, monkeypatch):
    source = tmp_path / "accepted_high_temperature_restart"
    source.write_bytes(b"analytical restart fixture, not a real wavefunction")
    captured = {}
    def complete(command, **kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(command, 0, "output fixture", "")
    monkeypatch.setattr(backend.subprocess, "run", complete)
    folder = tmp_path / "cooled"
    backend._execute(Atoms("He"), folder, "unused.exe", 0, 0, 1, ["--grad"], 1,
                     restart_file=source)
    assert "--restart" in captured["command"]
    assert (folder / "xtbrestart").read_bytes() == source.read_bytes()
    assert "sha256" in (folder / "launch.json").read_text()
