"""Full molecular RRHO thermochemistry with Grimme entropy interpolation.

Public energies are kcal/mol, entropies cal/(mol K), frequencies cm^-1.
The ``298`` field names follow the project schema; the actual state is ALWAYS
``temperature``. Enthalpy and ZPVE are retained from the supplied harmonic
calculation. Only vibrational entropy is replaced. ``sol`` means a concentration
standard state, and does not imply that a solvation free energy was computed.

References: Grimme, Chem. Eur. J. 2012, 18, 9955, doi:10.1002/chem.201200497;
https://gaussian.com/wp-content/uploads/dl/thermo.pdf;
https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Iterable, Literal

from pydantic import ConfigDict, Field
from pydantic.dataclasses import dataclass

R_J = 8.31446261815324
R_CAL = R_J / 4.184
H_PLANCK = 6.62607015e-34
K_B = 1.380649e-23
C_CM = 2.99792458e10
HARTREE_TO_KCAL = 627.5094740631
B_AV = 1.0e-44
_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?"
RotorInertia = Literal["grimme", "uncapped"]


@dataclass(frozen=True, kw_only=True, config=ConfigDict(allow_inf_nan=False))
class ThermochemistryResult:
    """Validated result, preserving full molecular entropy and provenance."""

    E_elec: float
    ZPVE: float = Field(ge=0)
    H_298: float
    S_harmonic: float = Field(ge=0)
    S_qRRHO: float
    G_298_qRRHO_sol: float
    temperature: float = Field(gt=0)
    pressure_pa: float = Field(gt=0)
    concentration_mol_l: float = Field(gt=0)
    cutoff_cm1: float = Field(gt=0)
    rotor_inertia: RotorInertia
    frequencies_cm1: tuple[float, ...]
    S_vib_harmonic: float
    S_vib_qRRHO: float
    delta_G_concentration: float
    G_qRRHO_gas: float
    source: str
    warnings: tuple[str, ...] = ()
    source_sha256: str | None = None
    program: str = "manual"
    job_index: int | None = None
    entropy_provenance: str = "User-supplied full molecular harmonic entropy"
    energy_unit: Literal["kcal/mol"] = "kcal/mol"
    entropy_unit: Literal["cal/(mol K)"] = "cal/(mol K)"
    frequency_unit: Literal["cm^-1"] = "cm^-1"
    temperature_unit: Literal["K"] = "K"
    # The source electronic energy may itself have used a continuum solvent.
    solvation_correction_added: bool = False


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return value


def _modes(frequencies: Iterable[float], *, transition_state: bool = False,
           exclude_zero_modes: bool = False) -> tuple[tuple[float, ...], list[str]]:
    values = tuple(float(v) for v in frequencies)
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError("Frequencies must be a nonempty sequence of finite cm^-1 values")
    negative = sum(v < 0 for v in values)
    if negative != int(transition_state):
        raise ValueError("A minimum must have no imaginary modes; transition_state=True "
                         "requires exactly one imaginary mode, which is excluded")
    zeros = values.count(0.0)
    if zeros and not exclude_zero_modes:
        raise ValueError("Zero frequencies require explicit exclude_zero_modes=True "
                         "for identified external translational/rotational modes")
    kept = tuple(v for v in values if v > 0)
    if not kept:
        raise ValueError("At least one positive vibrational frequency is required")
    warnings = []
    if negative:
        warnings.append("One imaginary transition-state mode excluded; no absolute-value conversion")
    if zeros:
        warnings.append(f"{zeros} identified zero external modes excluded")
    return kept, warnings


def _harmonic_mode(frequency: float, temperature: float) -> float:
    log_x = math.log(H_PLANCK * C_CM / K_B) + math.log(frequency) - math.log(temperature)
    if log_x < -9.0:
        x = math.exp(log_x)
        return R_CAL * (1.0 - log_x + x * x / 24.0)
    if log_x > math.log(700.0):
        return 0.0
    x = math.exp(log_x)
    return R_CAL * (x / math.expm1(x) - math.log(-math.expm1(-x)))


def harmonic_vibrational_entropy(frequencies_cm1: Iterable[float],
                                 temperature: float = 298.15) -> float:
    """Sum positive harmonic-mode entropies in cal/(mol K); reject zero/imaginary modes."""
    temperature = _positive(temperature, "temperature")
    modes, _ = _modes(frequencies_cm1)
    return math.fsum(_harmonic_mode(v, temperature) for v in modes)


def qrrho_vibrational_entropy(frequencies_cm1: Iterable[float], temperature: float = 298.15,
                              cutoff_cm1: float = 100.0,
                              rotor_inertia: RotorInertia = "grimme") -> float:
    """Grimme entropy in cal/(mol K), using mu_eff=mu*B_av/(mu+B_av).

    ``uncapped`` reproduces the directive's literal mu=h/(8*pi^2*nu), which
    diverges as nu approaches zero and is not the original reduced-inertia model.
    No frequency floor is applied. Computation uses logarithms for stability.
    """
    temperature = _positive(temperature, "temperature")
    cutoff_cm1 = _positive(cutoff_cm1, "cutoff_cm1")
    if rotor_inertia not in ("grimme", "uncapped"):
        raise ValueError("rotor_inertia must be 'grimme' or 'uncapped'")
    modes, _ = _modes(frequencies_cm1)
    terms = []
    for frequency in modes:
        log_mu = math.log(H_PLANCK / (8 * math.pi**2 * C_CM)) - math.log(frequency)
        if rotor_inertia == "grimme":
            z = math.log(B_AV) - log_mu
            log_mu = math.log(B_AV) - (max(z, 0.0) + math.log1p(math.exp(-abs(z))))
        log_q = 0.5 * (math.log(8 * math.pi**3 * K_B / H_PLANCK**2)
                       + log_mu + math.log(temperature))
        rotor = R_CAL * (0.5 + log_q)
        z = 4 * (math.log(frequency) - math.log(cutoff_cm1))
        e = math.exp(-abs(z))
        weight, rotor_weight = ((1 / (1 + e), e / (1 + e)) if z >= 0
                                else (e / (1 + e), 1 / (1 + e)))
        terms.append(weight * _harmonic_mode(frequency, temperature) + rotor_weight * rotor)
    return math.fsum(terms)


def concentration_correction(temperature: float = 298.15, pressure_pa: float = 101325.0,
                             concentration_mol_l: float = 1.0) -> float:
    """Ideal-gas p -> C standard-state shift RT ln(C*1000*RT/p), kcal/mol."""
    temperature = _positive(temperature, "temperature")
    pressure_pa = _positive(pressure_pa, "pressure_pa")
    concentration_mol_l = _positive(concentration_mol_l, "concentration_mol_l")
    log_ratio = (math.log(concentration_mol_l) + math.log(1000 * R_J)
                 + math.log(temperature) - math.log(pressure_pa))
    correction = (R_J / 4184) * temperature * log_ratio
    if not math.isfinite(correction):
        raise ValueError("Concentration correction exceeds finite floating-point range")
    return correction


def calculate_thermochemistry(
    frequencies_cm1: Iterable[float], *, E_elec: float, ZPVE: float,
    H_298: float, S_harmonic: float, temperature: float = 298.15,
    pressure_pa: float = 101325.0, concentration_mol_l: float = 1.0,
    cutoff_cm1: float = 100.0, rotor_inertia: RotorInertia = "grimme",
    transition_state: bool = False, exclude_zero_modes: bool = False,
    source: str = "manual", source_sha256: str | None = None,
    program: str = "manual", job_index: int | None = None,
    entropy_provenance: str = "User-supplied full molecular harmonic entropy",
    warnings: Iterable[str] = (),
) -> ThermochemistryResult:
    """Replace only the harmonic vibrational entropy in a complete RRHO result.

    E_elec, ZPVE, H_298: kcal/mol; S_harmonic: FULL molecular cal/(mol K).
    Inputs must refer to the same species, temperature, pressure and calculation.
    H already includes ZPVE: it is never added a second time. Solvation energy,
    conformational averaging and excited-state populations are not introduced.
    """
    temperature = _positive(temperature, "temperature")
    values = (E_elec, ZPVE, H_298, S_harmonic)
    if not all(math.isfinite(float(v)) for v in values) or ZPVE < 0 or S_harmonic < 0:
        raise ValueError("Energies/entropy must be finite; ZPVE and harmonic entropy nonnegative")
    modes, notes = _modes(frequencies_cm1, transition_state=transition_state,
                          exclude_zero_modes=exclude_zero_modes)
    s_harm = harmonic_vibrational_entropy(modes, temperature)
    if S_harmonic < s_harm - 0.1:
        raise ValueError("Full molecular harmonic entropy is smaller than its vibrational "
                         "component; check units, frequency scaling and entropy provenance")
    s_quasi = qrrho_vibrational_entropy(modes, temperature, cutoff_cm1, rotor_inertia)
    total_quasi = S_harmonic - s_harm + s_quasi
    shift = concentration_correction(temperature, pressure_pa, concentration_mol_l)
    gas = H_298 - temperature * total_quasi / 1000
    notes.extend(warnings)
    notes.append("Concentration standard-state correction only; no solvation free energy added")
    if rotor_inertia == "uncapped":
        notes.append("Uncapped literal-directive rotor inertia selected; diverges at zero frequency")
    return ThermochemistryResult(
        E_elec=E_elec, ZPVE=ZPVE, H_298=H_298, S_harmonic=S_harmonic,
        S_qRRHO=total_quasi, G_298_qRRHO_sol=gas + shift, temperature=temperature,
        pressure_pa=pressure_pa, concentration_mol_l=concentration_mol_l,
        cutoff_cm1=cutoff_cm1, rotor_inertia=rotor_inertia, frequencies_cm1=modes,
        S_vib_harmonic=s_harm, S_vib_qRRHO=s_quasi, delta_G_concentration=shift,
        G_qRRHO_gas=gas, source=source, warnings=tuple(notes),
        source_sha256=source_sha256, program=program, job_index=job_index,
        entropy_provenance=entropy_provenance,
    )


def _number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _required(pattern: str, text: str, label: str) -> float:
    matches = re.findall(pattern, text, re.I | re.M)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {label} in the selected analysis; found {len(matches)}")
    return _number(matches[0])


def _select_job(text: str, program: str, job_index: int) -> tuple[str, int]:
    termination = (r"^.*Normal termination of Gaussian[^\n]*" if program == "Gaussian"
                   else r"^.*ORCA TERMINATED NORMALLY[^\n]*")
    ends = list(re.finditer(termination, text, re.M | re.I))
    if not ends:
        raise ValueError(f"No normal termination found for {program}")
    chunks, start = [], 0
    for end in ends:
        chunks.append(text[start:end.end()])
        start = end.end()
    tail = text[start:]
    if job_index == -1 and re.search(
        r"SCF Done|Frequencies\s+--|Entering Gaussian|Error termination|"
        r"VIBRATIONAL FREQUENCIES|FINAL SINGLE POINT ENERGY|THERMOCHEMISTRY|"
        r"O\s+R\s+C\s+A|aborting the run|ORCA finished by error", tail, re.I
    ):
        raise ValueError("The final job is incomplete or failed; explicitly select an earlier job_index")
    try:
        block = chunks[job_index]
    except IndexError as exc:
        raise ValueError(f"job_index {job_index} is outside {len(chunks)} completed jobs") from exc
    index = job_index % len(chunks)
    boundaries = list(re.finditer(
        r"^\s*(?:--Link1--|Link1:\s+Proceeding to internal job step|\$new_job|"
        r"JOB NUMBER\s+\d+|Entering Gaussian System|\*\s+O\s+R\s+C\s+A\s+\*).*$",
        block, re.M | re.I,
    ))
    if boundaries:
        block = block[boundaries[-1].end():]
    if re.search(r"Error termination|ORCA finished by error|aborting the run", block, re.I):
        raise ValueError("Selected job contains an error termination; no cross-job recovery is attempted")
    return block, index


def _gaussian(block: str) -> dict:
    headers = list(re.finditer(r"^\s*Harmonic frequencies\b", block, re.M | re.I))
    if headers:
        block = block[headers[-1].start():]
    freq_lines = re.findall(r"^[ \t]*Frequencies[ \t]+--+[ \t]*([^\n]*)", block, re.M)
    if not freq_lines:
        raise ValueError("No Gaussian harmonic frequency block found")
    tokens = [token for line in freq_lines for token in line.split()]
    if any(not line.strip() for line in freq_lines) or any(re.fullmatch(_NUM, token) is None for token in tokens):
        raise ValueError("Malformed Gaussian frequency token; no damaged mode can be silently discarded")
    frequencies = [_number(token) for token in tokens]
    state = re.findall(rf"Temperature\s+({_NUM})\s+Kelvin\.\s+Pressure\s+({_NUM})\s+Atm", block, re.I)
    if len(state) != 1:
        raise ValueError("Gaussian analysis must contain one temperature/pressure record")
    zp = _required(rf"^\s*Zero-point correction=\s*({_NUM})", block, "ZPVE")
    ezp = _required(rf"^\s*Sum of electronic and zero-point Energies=\s*({_NUM})", block, "E+ZPVE")
    enthalpy = _required(rf"^\s*Sum of electronic and thermal Enthalpies=\s*({_NUM})", block, "enthalpy")
    # Anchor on the units header, rather than unrelated 'Total' lines elsewhere.
    table = re.search(r"E\s*\(Thermal\)\s+CV\s+S\s*\n\s*KCal/Mol\s+Cal/Mol-Kelvin\s+Cal/Mol-Kelvin(?P<table>[\s\S]+)", block, re.I)
    if not table:
        raise ValueError("Gaussian full harmonic entropy table is missing")
    entropy = _required(rf"^\s*Total\s+{_NUM}\s+{_NUM}\s+({_NUM})\s*$", table["table"], "total entropy")
    vibrational = re.findall(rf"^\s*Vibrational\s+{_NUM}\s+{_NUM}\s+({_NUM})\s*$", table["table"], re.M)
    if vibrational and any(v > 0 for v in frequencies):
        expected = harmonic_vibrational_entropy([v for v in frequencies if v > 0], _number(state[0][0]))
        if len(vibrational) != 1 or not math.isclose(_number(vibrational[0]), expected, abs_tol=0.03):
            raise ValueError("Gaussian printed vibrational entropy disagrees with parsed frequencies; "
                             "check scaling or incomplete/nonharmonic analysis")
    if re.search(r"Hindered Rotor|Anharmonic thermochemistry", block, re.I):
        raise ValueError("Hindered-rotor or anharmonic thermochemistry is not a harmonic RRHO baseline")
    return dict(frequencies_cm1=frequencies, E_elec=(ezp - zp) * HARTREE_TO_KCAL,
                ZPVE=zp * HARTREE_TO_KCAL, H_298=enthalpy * HARTREE_TO_KCAL,
                S_harmonic=entropy, temperature=_number(state[0][0]),
                pressure_pa=_number(state[0][1]) * 101325,
                entropy_provenance="Gaussian printed full molecular harmonic entropy table")


def _orca(block: str, transition_state: bool) -> dict:
    declared_cutoffs = re.findall(rf"\bCutOffFreq\s+({_NUM})", block, re.I)
    source_cutoff = _number(declared_cutoffs[-1]) if declared_cutoffs else 1.0
    headers = list(re.finditer(r"^\s*VIBRATIONAL FREQUENCIES\s*$", block, re.M))
    if not headers:
        raise ValueError("No ORCA vibrational frequency block found")
    block = block[headers[-1].end():]
    thermo = list(re.finditer(r"^\s*THERMOCHEMISTRY AT\b", block, re.M))
    if len(thermo) != 1:
        raise ValueError("ORCA parser requires one thermochemistry temperature per frequency analysis")
    freq_block, thermal = block[:thermo[0].start()], block[thermo[0].start():]
    freq_table = re.split(r"^\s*(?:NORMAL MODES|IR SPECTRUM|RAMAN SPECTRUM)\s*$", freq_block, maxsplit=1, flags=re.M)[0]
    candidate_rows = re.findall(r"^\s*\d+\s*:[^\n]*", freq_table, re.M)
    records = []
    for row in candidate_rows:
        match = re.fullmatch(rf"\s*(\d+)\s*:\s*({_NUM})\s*cm(?:\*\*|\^)?-1([^\n]*)", row, re.I)
        if match is None:
            raise ValueError("Malformed ORCA indexed frequency row; no damaged mode can be silently discarded")
        records.append(match.groups())
    if not records:
        raise ValueError("No ORCA indexed cm^-1 frequencies found")
    indexes = [int(record[0]) for record in records]
    if indexes != list(range(len(records))):
        raise ValueError("ORCA frequency indexes are incomplete, duplicated or out of order")
    frequencies = []
    for _, value, suffix in records:
        frequency = _number(value)
        if "imaginary" in suffix.lower() and frequency >= 0:
            raise ValueError("Imaginary-mode annotation conflicts with frequency sign")
        frequencies.append(frequency)
    zeros = [i for i, value in enumerate(frequencies) if value == 0]
    if zeros and (len(zeros) not in (5, 6) or zeros != list(range(len(zeros)))):
        raise ValueError("ORCA zero modes must be the initial five/six projected external modes")
    modes, _ = _modes(frequencies, transition_state=transition_state, exclude_zero_modes=True)
    if any(v <= source_cutoff for v in modes):
        raise ValueError("ORCA frequencies reach its thermochemistry CutOffFreq; rerun with a lower "
                         "source cutoff so printed enthalpy and all positive vibrational modes agree")
    included = [_number(v) for v in re.findall(rf"^\s*freq\.\s+({_NUM})\s+E\(vib\)", thermal, re.M)]
    if included and (len(included) != len(modes) or any(
        not math.isclose(a, b, rel_tol=0, abs_tol=0.02) for a, b in zip(included, modes)
    )):
        raise ValueError("ORCA thermal vibration list disagrees with the frequency block; "
                         "filtered thermochemistry cannot be silently combined with all modes")
    t = _required(rf"^\s*Temperature\s+\.\.\.\s*({_NUM})\s+K", thermal, "temperature")
    p = _required(rf"^\s*Pressure\s+\.\.\.\s*({_NUM})\s+atm", thermal, "pressure")
    def eh(label: str) -> float:
        return _required(rf"^\s*{label}\s+\.\.\.\s*({_NUM})\s+Eh", thermal, label)
    # ORCA can print Total enthalpy again in its Gibbs summary; isolate inner-energy/enthalpy first.
    before_entropy = re.split(r"^\s*ENTROPY\s*$", thermal, maxsplit=1, flags=re.M)[0]
    h = _required(rf"^\s*Total Enthalpy\s+\.\.\.\s*({_NUM})\s+Eh", before_entropy, "enthalpy")
    # These components are unaffected by ORCA's default qRRHO, preventing double damping.
    try:
        nonvib_ts = sum(eh(label + " entropy") for label in ("Electronic", "Rotational", "Translational"))
    except ValueError as exc:
        raise ValueError("ORCA requires electronic, rotational and translational entropy components "
                         "to reconstruct harmonic entropy safely; its total may already use qRRHO") from exc
    entropy = nonvib_ts * HARTREE_TO_KCAL * 1000 / t + harmonic_vibrational_entropy(modes, t)
    return dict(frequencies_cm1=frequencies, E_elec=eh("Electronic energy") * HARTREE_TO_KCAL,
                ZPVE=eh("Zero point energy") * HARTREE_TO_KCAL, H_298=h * HARTREE_TO_KCAL,
                S_harmonic=entropy, temperature=t, pressure_pa=p * 101325,
                exclude_zero_modes=True,
                entropy_provenance="ORCA nonvibrational entropy components plus recomputed harmonic vibration",
                warnings=("ORCA printed vibrational/total entropy deliberately bypassed to avoid double qRRHO correction",))


def parse_thermochemistry(path: str | Path, *, job_index: int = -1,
                         temperature: float | None = None, transition_state: bool = False,
                         concentration_mol_l: float = 1.0, cutoff_cm1: float = 100.0,
                         rotor_inertia: RotorInertia = "grimme") -> ThermochemistryResult:
    """Parse a normally terminated Gaussian 16 .log or ORCA .out directly.

    Each normal termination defines a completed job (zero-based ``job_index``).
    Default: last completed job, rejecting a subsequently failed/incomplete run.
    Within that job the last frequency analysis is used. Records from different
    jobs are never combined. A temperature override must equal the parsed state:
    changing it requires recomputation of all RRHO contributions upstream.
    """
    path = Path(path)
    raw = path.read_bytes()
    text = raw.decode("utf-8-sig", errors="replace")
    gaussian = bool(re.search(r"Normal termination of Gaussian|Frequencies\s+--", text))
    orca = bool(re.search(r"ORCA TERMINATED NORMALLY|VIBRATIONAL FREQUENCIES", text))
    if gaussian == orca:
        raise ValueError("Cannot unambiguously identify Gaussian or ORCA output")
    program = "Gaussian" if gaussian else "ORCA"
    block, index = _select_job(text, program, job_index)
    values = _gaussian(block) if gaussian else _orca(block, transition_state)
    if temperature is not None and not math.isclose(
        _positive(temperature, "temperature"), values["temperature"], abs_tol=1e-6, rel_tol=0
    ):
        raise ValueError("Requested temperature differs from parsed full RRHO state; rerun thermochemistry upstream")
    return calculate_thermochemistry(
        **values, transition_state=transition_state, concentration_mol_l=concentration_mol_l,
        cutoff_cm1=cutoff_cm1, rotor_inertia=rotor_inertia, source=str(path.resolve()),
        source_sha256=hashlib.sha256(raw).hexdigest(), program=program, job_index=index,
    )
