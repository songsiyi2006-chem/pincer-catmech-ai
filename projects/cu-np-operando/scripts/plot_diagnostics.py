"""Plot actual compute status and a clearly labeled numerical convergence audit."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / 'results/figures'
    out.mkdir(exist_ok=True, parents=True)
    audit = json.loads((root / 'results/neutral_transport_audit.json').read_text())
    grid = audit['grid_verification']['grids']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    n = [v['cells'] for v in grid]
    err = [v['absolute_error_mol_m3'] for v in grid]
    ax.loglog(n, err, 'o-', color='#176B87', label='Measured grid error vs analytic solution')
    ax.loglog(n, [err[0] * n[0] / x for x in n], '--', color='#B86E30', label='First-order reference')
    ax.set(xlabel='Axial cells', ylabel='Outlet concentration error (mol / m3)',
           title='Numerical verification only — neutral synthetic fixture')
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / 'neutral_grid_convergence.png', dpi=180)
    plt.close(fig)
    summary = json.loads((root / 'data/dft_pilot/summary.json').read_text())
    runs = [r for r in summary['runs'] if r.get('elapsed_seconds') is not None]
    if runs:
        fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
        labels = [r['run_id'] for r in runs]
        colors = ['#176B87' if r['status'] == 'completed' else '#B86E30' for r in runs]
        bars = ax.barh(labels, [r['elapsed_seconds'] for r in runs], color=colors)
        for bar, r in zip(bars, runs):
            ax.text(bar.get_width()+3, bar.get_y()+bar.get_height()/2, r['status'], va='center', fontsize=8)
        ax.set_xlim(0, max(r['elapsed_seconds'] for r in runs)*1.35)
        ax.set(xlabel='Observed wall time (s)', title='Actual bounded molecular DFT jobs — not electrode calculations')
        fig.savefig(out / 'native_job_status.png', dpi=180)
        plt.close(fig)


if __name__ == '__main__':
    main()
