"""Plot actual archived thermochemical sensitivity; ranges are not confidence intervals."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]


def main():
    with (REPO/'data/phase4/solvation_sensitivity_001/sensitivity.csv').open(encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    rows = [{k:(v if k=='name' else float(v)) for k,v in r.items()} for r in rows]
    colors = ['#0072B2','#D55E00','#009E73']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,axes = plt.subplots(1,2,figsize=(9.3,4.2),layout='constrained')
    for n,color in zip((1,2,3),colors):
        cut = sorted([r for r in rows if r['n_alcohol']==n and r['temperature_K']==383.15
                      and r['alcohol_activity']==1],key=lambda r:r['cutoff_cm1'])
        act = sorted([r for r in rows if r['n_alcohol']==n and r['temperature_K']==383.15
                      and r['cutoff_cm1']==100],key=lambda r:r['alcohol_activity'])
        axes[0].plot([r['cutoff_cm1'] for r in cut],[r['delta_G_standard_1M_kcal_mol'] for r in cut],
                     'o-',color=color,label=f'n = {n}')
        axes[1].plot([r['alcohol_activity'] for r in act],[r['effective_association_kcal_mol'] for r in act],
                     'o-',color=color,label=f'n = {n}')
    axes[0].set(xlabel=r'qRRHO cutoff (cm$^{-1}$)',ylabel='Association free energy (kcal/mol)',
                title='A  Standard association, 1 M',xticks=[50,100,150])
    axes[1].set(xlabel='Neutral tBuOH activity (relative to 1 M)',ylabel='Conditional free energy (kcal/mol)',
                title='B  Activity sensitivity, cutoff 100',xscale='log')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.legend(frameon=False)
        ax.grid(alpha=.15)
    fig.suptitle('Selected xTB/ALPB minima at 383.15 K; no bulk population or barrier inference',fontsize=11)
    folder=REPO/'examples/plots'
    fig.savefig(folder/'phase4_solvation_sensitivity.png',dpi=240)
    fig.savefig(folder/'phase4_solvation_sensitivity.svg')
    svg = folder/'phase4_solvation_sensitivity.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',
                   encoding='utf-8', newline='\n')
    plt.close(fig)


if __name__=='__main__':
    main()
