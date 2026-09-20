"""Auditable *neutral* three-state software verification kernel.

This is not a Cu electrosynthesis reaction mechanism. Energies in examples are
deliberately synthetic, and the module does not calculate an electrode potential,
current, electrochemical free energy, or catalyst ranking.
"""

from .core import (
    CycleEnergies, RateConstants, SurfaceState, FilmState, PlugFlowState,
    ElectricalAccounting, InputEvidence, rate_constants, steady_surface,
    film_flux, plug_flow, plug_flow_adaptive, cell_energy, faradaic_efficiency,
    target_ranking_gate, K_B_EV_K, K_B_J_K, PLANCK_J_S, FARADAY_C_MOL,
)

__all__ = [name for name in globals() if not name.startswith("_")]
