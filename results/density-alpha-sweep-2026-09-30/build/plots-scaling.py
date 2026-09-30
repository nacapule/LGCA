# Divergence test: late cluster density (plateau) vs total particle number at FIXED average density 2.4 (input .4),
# box sizes 90/120/180/240, centre + neighbours. One line per alpha; bars = seed range. Also time curves per alpha.
# A (alpha, L) point uses every seed that reached T (at least 2). Late = per-run mean over steps 0.75T+100..T.
# Usage: python3 build/plots-scaling.py [T]
import json, sys, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
out = Path(__file__).resolve().parent.parent
figdir = out / 'figures'
T = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
sizes = [90, 120, 180, 240]; alphas = [.8, .82, .84, .85, .86, .88, .9, 1]; seeds = [12345, 777, 424242]
lo = int(T * .75) + 100
def load(L, a, s):
    p = out / 'data' / f'both-L{L}-a{a:g}-d0.4-seed{s}' / 'run.json'
    if not p.exists(): return None
    r = json.loads(p.read_text())
    if r['snaps'][-1]['t'] < T: return None
    ts = np.array([m['t'] for m in r['series'] if m['t'] <= T]); v = np.array([m['clusterDensity'] for m in r['series'] if m['t'] <= T])
    return ts, v, r['series'][0]['N']
rows = []; curves = {}
for a in alphas:
    for L in sizes:
        rs = [x for x in (load(L, a, s) for s in seeds) if x]
        if len(rs) < 2: continue
        late = [x[1][x[0] >= lo].mean() for x in rs]
        rows.append({'alpha': a, 'L': L, 'N': float(np.mean([x[2] for x in rs])), 'seeds': len(rs), 'T': T,
                     'plateau': float(np.mean(late)), 'seed_min': float(min(late)), 'seed_max': float(max(late))})
        curves[(a, L)] = (rs[0][0], np.mean([x[1] for x in rs], axis=0))
with open(out / f'scaling-t{T}.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': '#fcfcfa',
                     'axes.facecolor': '#fcfcfa', 'savefig.facecolor': '#fcfcfa', 'text.color': '#19252d', 'axes.labelcolor': '#19252d'})
kfmt = FuncFormatter(lambda v, _: f'{v:g}' if v < 1000 else f'{v / 1000:g}k')
present = [a for a in alphas if any(r['alpha'] == a for r in rows)]
acm = plt.get_cmap('plasma'); acol = {a: acm(0.85 * i / max(1, len(present) - 1)) for i, a in enumerate(present)}

# Figure A: plateau vs N, the divergence test.
fig, ax = plt.subplots(figsize=(7.6, 5), layout='constrained')
for a in present:
    rr = sorted([r for r in rows if r['alpha'] == a], key=lambda r: r['N'])
    x = [r['N'] for r in rr]; y = [r['plateau'] for r in rr]
    ax.errorbar(x, y, yerr=[[r['plateau'] - r['seed_min'] for r in rr], [r['seed_max'] - r['plateau'] for r in rr]],
                fmt='-o', ms=4, lw=1.6, capsize=3, color=acol[a], label=f'α {a:g}')
ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(alpha=.18, which='both')
ax.xaxis.set_major_formatter(kfmt); ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_major_formatter(kfmt)
Ns = sorted({round(r['N'], -2) for r in rows})
ax.set_xticks([r['N'] for r in rows if r['alpha'] == present[-1]] or Ns)
ax.set_xticklabels([f"{r['N'] / 1000:.0f}k\n(L {r['L']})" for r in rows if r['alpha'] == present[-1]] or [f'{n / 1000:.0f}k' for n in Ns])
ax.set_xlabel('Total particles N (box size L), average density fixed at 2.4 particles/site')
ax.set_ylabel('Plateau cluster density (particles/site)')
ax.set_title(f'Flat = bounded, rising = diverging · steps {lo:,}–{T:,} · bar = seed range',
             loc='left', fontsize=9.5, color='#5a676c')
ax.legend(frameon=False, fontsize=8, ncol=2, loc='upper left')
for ext in ['png', 'pdf']: fig.savefig(figdir / f'scaling-plateau-vs-N-t{T}.{ext}', dpi=150)
plt.close(fig)

# Figure B: time curves per alpha, one line per box size.
scol = {90: '#c6dbef', 120: '#6baed6', 180: '#2171b5', 240: '#08306b'}
n = len(present); ncol = 4; nrow = -(-n // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.7 * nrow + .5), layout='constrained', squeeze=False, sharey=True)
for ax, a in zip(axes.flat, present):
    for L in sizes:
        if (a, L) in curves:
            ts, v = curves[(a, L)]; ax.plot(ts[1:], v[1:], color=scol[L], lw=1.1, label=f'L {L}')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(25, T); ax.grid(alpha=.18)
    ax.xaxis.set_major_formatter(kfmt); ax.yaxis.set_major_formatter(kfmt); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(f'α = {a:g}', loc='left', fontsize=10, fontweight='bold')
for ax in axes.flat[n:]: ax.axis('off')
axes[0, 0].legend(frameon=False, fontsize=8, loc='lower right')
for ax in axes[:, 0]: ax.set_ylabel('Cluster density')
for ax in axes[-1, :]: ax.set_xlabel('Simulation step')
fig.suptitle(f'Cluster density vs time at fixed average density 2.4, by box size · centre + neighbours, sens 6 · mean of 2–3 seeds',
             x=.01, ha='left', fontsize=10, color='#5a676c')
for ext in ['png', 'pdf']: fig.savefig(figdir / f'scaling-time-t{T}.{ext}', dpi=150)
plt.close(fig)
for r in rows: print(f"α {r['alpha']:<5} L{r['L']:<4} N {r['N']:8.0f} seeds {r['seeds']}  plateau {r['plateau']:7.1f}  ({r['seed_min']:.1f}–{r['seed_max']:.1f})")
