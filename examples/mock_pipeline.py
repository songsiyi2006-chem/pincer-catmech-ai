"""Runnable synthetic Ru-PNP-like demonstration, not an optimized molecule.

The eight-atom toy scaffold has 18 assigned positive modes (3N-6). Energies
and the 65 cal/(mol K) non-vibrational entropy are explicitly artificial.
"""

from dataclasses import asdict
import json

from ase import Atoms

from pincer_catmech.features.steric import bite_angle, buried_volume
from pincer_catmech.kinetics.free_energy import (
    calculate_thermochemistry, harmonic_vibrational_entropy,
)

FREQUENCIES = (
    25.4, 68.2, 120.5, 160.0, 210.0, 300.0, 420.0, 550.0, 700.0,
    850.0, 1000.0, 1100.0, 1250.0, 1400.0, 1500.0, 1650.0, 3000.0, 3100.0,
)


def run_demo():
    """Return JSON-serializable numerical results and explicit data provenance."""
    atoms = Atoms(
        symbols=["Ru", "P", "P", "N", "C", "C", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0], [2.3, 0.0, 0.0], [-2.3, 0.0, 0.0],
            [0.0, 2.1, 0.0], [1.6, 1.5, 0.0], [-1.6, 1.5, 0.0],
            [1.7, 1.5, 1.1], [-1.7, 1.5, 1.1],
        ],
    )
    result = calculate_thermochemistry(
        FREQUENCIES, E_elec=-1000000.0, ZPVE=25.0, H_298=-999970.0,
        S_harmonic=65.0 + harmonic_vibrational_entropy(FREQUENCIES),
        source="Synthetic Ru-PNP-like unit fixture; not quantum chemistry output",
    )
    profile = buried_volume(atoms, 0, ligand_indices=range(1, len(atoms)), seed=2026)
    return {
        "evidence_level": "synthetic software demonstration",
        "charge_and_spin": "not assigned; no electronic structure calculation",
        "thermochemistry": asdict(result),
        "bite_angle_P_Ru_P_degrees": bite_angle(atoms, 0, 1, 2),
        "bite_angle_P_Ru_N_degrees": bite_angle(atoms, 0, 1, 3),
        "steric": asdict(profile),
    }


if __name__ == "__main__":
    print(json.dumps(run_demo(), ensure_ascii=False, indent=2, allow_nan=False))
