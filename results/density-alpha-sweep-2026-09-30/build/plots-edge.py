# Edge profiles figure from edge-profiles-t{T}.json (build/edge.mjs). Negative distance = inside the largest cluster.
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out = Path(__file__).resolve().parent.parent; T = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
rows = [r for r in json.loads((out / f'edge-profiles-t{T}.json').read_text()) if r['seeds'] == 5 and r['alpha'] in (.8, .9, 1)]
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': '#fcfcfa',
                     'axes.facecolor': '#fcfcfa', 'savefig.facecolor': '#fcfcfa', 'text.color': '#19252d', 'axes.labelcolor': '#19252d'})
dens = sorted({r['dens'] for r in rows}); alphas = sorted({r['alpha'] for r in rows})
cols = {a: c for a, c in zip(alphas, ['#d7191c', '#2c7bb6', '#084081'])}
fig, axes = plt.subplots(2, len(dens), figsize=(4.2 * len(dens), 6.2), layout='constrained', sharex=True)
for j, d in enumerate(dens):
    for r in [r for r in rows if r['dens'] == d]:
        p = r['profile']; x = [b['distance'] for b in p]
        axes[0, j].plot(x, [b['density'] / r['mean'] for b in p], '-o', ms=2.5, color=cols[r['alpha']], label=f"α {r['alpha']:g}")
        q = [b for b in p if b['localPolar'] is not None]
        axes[1, j].plot([b['distance'] for b in q], [b['localPolar'] for b in q], '-o', ms=2.5, color=cols[r['alpha']])
    axes[0, j].set_yscale('log'); axes[0, j].axhline(1, color='#999', lw=.8, ls=':')
    axes[0, j].set_title(f'average density {6 * d:g} particles/site', loc='left', fontsize=10, fontweight='bold')
    for ax in axes[:, j]: ax.axvline(0, color='#999', lw=.8, ls='--'); ax.grid(alpha=.18)
    axes[1, j].set_xlabel('Hops from cluster edge  (← inside · outside →)'); axes[1, j].set_ylim(0, 1.02)
axes[0, 0].set_ylabel('Density ÷ average density'); axes[1, 0].set_ylabel('Local polar order (occupied sites)')
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle(f'Density and local order vs distance to the largest cluster\'s edge · t = {T:,} · 90×90, centre + neighbours · mean of 5 seeds',
             x=.01, ha='left', fontsize=10, color='#5a676c')
for ext in ['png', 'pdf']: fig.savefig(out / 'figures' / f'edge-profiles-t{T}.{ext}', dpi=150)
