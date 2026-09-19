"""Small evidence-bound scientific utilities; eV per event unless stated."""
import math
import re
import hashlib

HARTREE_EV = 27.211386245981
EV_KJMOL = 96.48533212331002
KB_EV = 8.617333262145e-5
H_EV_S = 4.135667696e-15
F_C_MOL = 96485.33212331002


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_xtb(text, task="sp"):
    """Reject energy-looking output without a completed native calculation."""
    values = re.findall(r"TOTAL ENERGY\s+([-+\d.Ee]+)\s+Eh", text)
    normal = "normal termination of xtb" in text.lower()
    failed = any(s in text.lower() for s in ["abnormal termination", "scc is not converged", "failed to converge"])
    optimized = "GEOMETRY OPTIMIZATION CONVERGED" in text
    if not values or not normal or failed or (task == "opt" and not optimized):
        raise ValueError("Incomplete, failed, or nonconverged xTB output")
    energy = float(values[-1])
    if not math.isfinite(energy):
        raise ValueError("nonfinite energy")
    return {"energy_hartree": energy, "native_normal_termination": normal,
            "geometry_optimized": optimized, "task": task}


def balanced(stoichiometry, species):
    """Products positive, reactants negative; explicit e- and H+ allowed."""
    totals = {}
    for name, coefficient in stoichiometry.items():
        for key, count in species[name].items():
            totals[key] = totals.get(key, 0) + coefficient * count
    return all(abs(x) < 1e-12 for x in totals.values()), totals


def convert_potential(value, source, target):
    """Calibration offsets against a shared reference; never infer nonaqueous pH."""
    required = ("solvent", "temperature_K", "common_reference", "offset_V", "source")
    if any(k not in x for x in (source, target) for k in required):
        raise ValueError("Explicit referenced calibration required")
    for key in ("solvent", "temperature_K", "common_reference"):
        if source[key] != target[key]:
            raise ValueError("Incompatible reference calibration")
    return value + source["offset_V"] - target["offset_V"]


def reservoir_free_energy(delta_g_eV, delta_n, chemical_potentials_eV):
    """Grand-state comparison: Delta Omega = Delta G - sum(Delta n_i mu_i)."""
    if set(delta_n) - set(chemical_potentials_eV):
        raise ValueError("Missing chemical reservoir")
    return delta_g_eV - sum(n * chemical_potentials_eV[k] for k, n in delta_n.items())


def rate_pair(g_reactant, g_product, g_ts, temperature_K=298.15):
    if temperature_K <= 0 or g_ts < max(g_reactant, g_product):
        raise ValueError("Positive temperature and a shared barrier above both states required")
    kt = KB_EV * temperature_K
    prefactor = kt / H_EV_S
    return (prefactor * math.exp(-(g_ts-g_reactant)/kt),
            prefactor * math.exp(-(g_ts-g_product)/kt))


def site_steady_state(rates):
    """Three sites: * <-> S* <-> P* <-> *. Rates include explicit activities.

    Input order adsorption +/-; conversion +/-; release +/-. This is an
    intrinsic, single-site, ideal Markov network, not an electrode current.
    """
    import numpy as np
    if len(rates) != 6 or min(rates) <= 0 or not np.isfinite(rates).all():
        raise ValueError("Six finite positive rates required")
    a,b,c,d,e,f = rates
    q = np.array([[-a-f,b,e], [a,-b-c,d], [f,c,-d-e]], dtype=float)
    rhs = np.array([0.,0.,1.])
    matrix = q.copy(); matrix[-1,:] = 1
    theta = np.linalg.solve(matrix,rhs)
    if min(theta) < -1e-10 or max(abs(q @ theta)) > 1e-8 * max(rates):
        raise ValueError("Invalid site balance")
    return theta, a*theta[0]-b*theta[1]


def grouped_split(records, holdout_groups):
    """Whole predeclared groups are reserved; repeated calculations stay together."""
    train, holdout = [], []
    group_by_identity = {}
    for row in records:
        old = group_by_identity.setdefault(row["identity"], row["group"])
        if old != row["group"]:
            raise ValueError("Same identity assigned to different groups")
        (holdout if row["group"] in holdout_groups else train).append(row)
    return train, holdout


def electrical_kwh_kg(full_cell_V, electrons, faradaic_efficiency, mw_g_mol, recovery=1.):
    """FE is moles of target formed per theoretical electron equivalent.

    Recovery is isolated acceptable target / target formed, not conversion.
    Auxiliary power and solvent/separation energy are excluded.
    """
    if full_cell_V <= 0 or electrons <= 0 or mw_g_mol <= 0:
        raise ValueError("Positive physical inputs required")
    if not (0 < faradaic_efficiency <= 1 and 0 < recovery <= 1):
        raise ValueError("FE and recovery must be fractions in (0,1]")
    return full_cell_V*electrons*F_C_MOL/(3600*mw_g_mol*faradaic_efficiency*recovery)
