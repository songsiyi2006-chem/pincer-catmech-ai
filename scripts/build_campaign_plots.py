"""Render measured structural/endpoint data and explicit unavailable kinetics.

No interpolated barrier, invented TOF, or zero-filled missing result is drawn.
SVG and PDF are vector outputs; PNG companions are used in the monographs.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "examples/plots"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save(fig, stem):
    for suffix in ("svg", "pdf", "png"):
        path = OUTPUT / f"{stem}.{suffix}"
        if suffix == "svg":
            with path.open("w", encoding="utf-8", newline="\n") as stream:
                fig.savefig(stream, format="svg", dpi=200, bbox_inches="tight", facecolor="white")
            path.write_text("\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n",
                            encoding="utf-8", newline="\n")
        else:
            fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    descriptors = REPO / "data/datasets/catalyst_descriptors.csv"
    reactions = REPO / "data/datasets/thermochemistry/reaction_free_energies.csv"
    conditions = REPO / "data/campaign/conditions.json"
    structure = rows(descriptors)
    thermal = rows(reactions)
    settings = json.loads(conditions.read_text())
    metals = ["Ru", "Mn", "Fe", "Co"]
    backbones = ["macho_pnp", "pyridine_pnn", "bipyridine_pnnoh"]
    columns = [(backbone, substituent) for backbone in backbones for substituent in ("Ph", "iPr")]
    values = np.full((4, 6), np.nan)
    for row in structure:
        if row.get("buried_volume_percent"):
            values[metals.index(row["metal"]), columns.index((row["backbone"], row["substituent"]))] = float(row["buried_volume_percent"])
    fig, ax = plt.subplots(figsize=(10, 4.8), layout="constrained")
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e5e7eb")
    im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=40, vmax=80, aspect="auto")
    ax.set_xticks(range(6), [f"{name}\n{r}" for name in ("PNP", "pyridine PNN", "bipyridine PNN(O)") for r in ("Ph", "iPr")])
    ax.set_yticks(range(4), ["Ru(II)", "Mn(I)", "Fe(II)", "Co(I)"])
    ax.set_title("Buried volume of identity-preserving active structures", pad=15, fontweight="bold")
    for i in range(4):
        for j in range(6):
            v = values[i, j]
            ax.text(j, i, f"{v:.2f}%" if np.isfinite(v) else "identity\nnot retained",
                    ha="center", va="center", fontsize=10,
                    color="white" if np.isfinite(v) and v < 65 else "#172b4d")
            if not np.isfinite(v):
                ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, fill=False, hatch="///", edgecolor="#b8bec7", linewidth=0))
    fig.colorbar(im, ax=ax, label="Ligand buried volume (%)", shrink=.85)
    fig.supxlabel("Best found electronic-energy conformers; 3.5 Å sphere, 1.17 × Bondi radii, H excluded.\n50,000 Monte Carlo points; Wilson sampling intervals are provided in the source CSV.", fontsize=9)
    save(fig, "steric_buried_volume_map")

    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True, layout="constrained")
    source_rows = []
    for ax, substituent in zip(axes, ("Ph", "iPr")):
        ax.axvspan(.40, .60, color="#eef0f3", zorder=0)
        for metal, color, offset in (("Ru", "#a33b63", -.06), ("Mn", "#127c80", .06)):
            identifier = f"{metal}_macho_pnp_{substituent}"
            selected = [row for row in thermal if row["catalyst_id"] == identifier
                        and row["reaction"] == "alcohol_dehydrogenation"
                        and abs(float(row["temperature_K"]) - 383.15) < 1e-8]
            if len(selected) != 1:
                continue
            row = selected[0]
            source_rows.append(row)
            delta = float(row["delta_G_kcal_mol"])
            ax.plot([.03+offset, .23+offset], [0, 0], lw=2.5, color=color)
            ax.plot([.77+offset, .97+offset], [delta, delta], lw=2.5, color=color,
                    label=f"{metal}: ΔG = {delta:+.2f} kcal/mol")
            ax.scatter([.13+offset, .87+offset], [0, delta], s=32, color=color)
        ax.set_xlim(-.09, 1.11)
        ax.set_xticks([.13, .87], ["Cat + alcohol", "Cat–H₂ + aldehyde"])
        ax.set_title(f"MACHO-type PNP, R = {substituent}")
        ax.text(.50, .95, "TS free energy\nunavailable", transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="#545f70")
        ax.legend(loc="lower left", fontsize=8, frameon=False)
        ax.grid(axis="y", alpha=.18)
    axes[0].set_ylabel("Reaction ΔG relative to separated reactants (kcal/mol)")
    fig.suptitle("Certified endpoint free energies at 383.15 K", fontweight="bold")
    fig.supxlabel("GFN2-xTB / ALPB toluene + qRRHO, 1 M. Each endpoint uses its own certified molecular minimum.\nThese separated-species thermodynamic levels are not a transition-state profile or an activation barrier.", fontsize=9)
    save(fig, "pes_comparison_ru_vs_mn")

    fig, ax = plt.subplots(figsize=(9, 5.6), layout="constrained")
    temperatures = settings["temperature_grid_K"]
    loading = np.array(settings["catalyst_loading_fraction_grid"]) * 100
    ax.set_xlim(min(temperatures)-5, max(temperatures)+5)
    ax.set_ylim(min(loading)-.05, max(loading)+.05)
    ax.add_patch(Rectangle((min(temperatures)-5, min(loading)-.05), 110, 2,
                          facecolor="#eff2f6", edgecolor="#d1d7df", hatch="///", linewidth=.5))
    ax.text(390, 1.14, "TOF NOT COMPUTED", ha="center", va="center", fontsize=20, color="#233f55", weight="bold")
    ax.text(390, .82, "Complete certified transition-state network and\nfree-base speciation are unavailable.",
            ha="center", va="center", fontsize=11, color="#334b5c", bbox={"facecolor":"#eff2f6", "edgecolor":"none", "pad":7})
    ax.set_xticks(temperatures)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Catalyst loading (mol%)")
    ax.set_title("Requested kinetic grid: result availability", fontweight="bold", pad=13)
    fig.supxlabel("Total tBuOK: 0.01–0.20 equiv; baseline 0.05 equiv. No numerical TOF values or color scale exist.\nHatching denotes unavailable calculations, not zero turnover. This is an availability plot, not a kinetic heatmap.", fontsize=9)
    save(fig, "microkinetic_tof_heatmap")
    manifest = {"sources": {str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in (descriptors, reactions, conditions)},
                "steric_finite_cells": int(np.isfinite(values).sum()), "steric_missing_cells": int(np.isnan(values).sum()),
                "endpoint_rows": source_rows, "temperature_K": 383.15,
                "pes_scope": "Certified separate-species endpoint thermodynamics; no TS curve",
                "heatmap_scope": "availability_plot_not_TOF_result", "numerical_TOF_values": 0,
                "missing_values_are_zero": False}
    (OUTPUT / "PLOT_PROVENANCE.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8", newline="\n")
    print(json.dumps({key: manifest[key] for key in ("steric_finite_cells", "steric_missing_cells", "heatmap_scope")}))


if __name__ == "__main__":
    main()
