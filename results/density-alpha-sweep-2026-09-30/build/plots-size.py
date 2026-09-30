# Box-size scaling with the centre included: density input .4 (2.4 particles/site), L = 90/120/180.
# L90 runs come from the main sweep (same settings). Uses runs that reached T; seeds common to all sizes of an alpha.
# Usage: python3 build/plots-size.py [T]
import json, sys, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
out = Path(__file__).resolve().parent.parent
figdir = out / 'figures'; figdir.mkdir(exist_ok=True)
T = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
sizes = [90, 120, 180]; alphas = [.8, .85, .9, 1]; seeds = [12345, 777, 424242]; lo = int(T * .75) + 100
def load(L, a, s):
    p = out / 'data' / f'both-L{L}-a{a:g}-d0.4-seed{s}' / 'run.json'
    if not p.exists(): return None
    r = json.loads(p.read_text())
    if r['snaps'][-1]['t'] < T: return None
    ts = np.array([m['t'] for m in r['series'] if m['t'] <= T]); v = np.array([m['clusterDensity'] for m in r['series'] if m['t'] <= T])
    return ts, v, r['series'][0]['N'] / L ** 2
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': '#fcfcfa',
                     'axes.facecolor': '#fcfcfa', 'savefig.facecolor': '#fcfcfa', 'text.color': '#19252d', 'axes.labelcolor': '#19252d'})
kfmt = FuncFormatter(lambda v, _: f'{v:g}' if v < 1000 else f'{v / 1000:g}k')
colors = {90: '#9ecae1', 120: '#4292c6', 180: '#08306b'}
rows = []; present = []
alphas = [a for a in alphas if [s for s in seeds if all(load(L, a, s) for L in sizes)]]
fig, axes = plt.subplots(1, len(alphas), figsize=(4 * len(alphas), 3.6), layout='constrained', sharey=True)
for ax, a in zip(axes, alphas):
    common = [s for s in seeds if all(load(L, a, s) for L in sizes)]
    if not common: ax.set_title(f'α = {a:g} · pending', loc='left'); continue
    present.append(a)
    for L in sizes:
        runs = [load(L, a, s) for s in common]; ts = runs[0][0]; v = np.mean([r[1] for r in runs], axis=0)
        late = [r[1][r[0] >= lo].mean() for r in runs]
        rows.append({'alpha': a, 'L': L, 'seeds': len(common), 'T': T, 'mean_density': float(np.mean([r[2] for r in runs])),
                     'late_cluster_density': float(np.mean(late)), 'seed_min': float(min(late)), 'seed_max': float(max(late))})
        ax.plot(ts[1:], v[1:], color=colors[L], lw=1.2, label=f'{L}×{L} · N ≈ {runs[0][2] * L * L:,.0f}')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(25, T); ax.grid(alpha=.18)
    ax.xaxis.set_major_formatter(kfmt); ax.yaxis.set_major_formatter(kfmt); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(f'α = {a:g} · {len(common)} seeds', loc='left', fontsize=10, fontweight='bold'); ax.set_xlabel('Simulation step')
axes[0].set_ylabel('Cluster density (particles/site)'); axes[0].legend(frameon=False, fontsize=8, loc='lower right')
fig.suptitle(f'Box size at fixed average density 2.4 particles/site · centre + neighbours, sens 6 · to {T:,} steps',
             x=.01, ha='left', fontsize=10, color='#5a676c')
for ext in ['png', 'pdf']: fig.savefig(figdir / f'size-time-t{T}.{ext}', dpi=150)
plt.close(fig)
if rows:
    with open(out / f'size-summary-t{T}.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(f"alpha {r['alpha']:<5} L{r['L']:<4} seeds {r['seeds']}  late {r['late_cluster_density']:8.1f}  (seeds {r['seed_min']:.1f}-{r['seed_max']:.1f})")
