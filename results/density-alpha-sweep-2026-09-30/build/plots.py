# Figures + summary table for the density x alpha sweep. Works on partial data: for horizon T it uses every run
# that has reached T (series cut at T). Late level = per-run mean over steps 0.75T+100..T, then mean over seeds.
# Usage: python3 build/plots.py [T]   (default: largest horizon reached by every seed-12345 job)
import json, sys, math, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from matplotlib.colors import LogNorm, TwoSlopeNorm
out = Path(__file__).resolve().parent.parent
figdir = out / 'figures'; figdir.mkdir(exist_ok=True)
A = 90 * 90
runs = []
for p in sorted((out / 'data').glob('both-L90-*/run.json')):
    r = json.loads(p.read_text())
    runs.append({'alpha': r['job']['alpha'], 'dens': r['job']['dens'], 'seed': r['job']['seed'],
                 'reached': r['snaps'][-1]['t'], 'series': r['series']})
if len(sys.argv) > 1:
    T = int(sys.argv[1])
else:
    T = min(r['reached'] for r in runs if r['seed'] == 12345)
use = [r for r in runs if r['reached'] >= T]
if len(sys.argv) > 2:  # optional seed filter, e.g. 12345,777
    keep = {int(x) for x in sys.argv[2].split(',')}; use = [r for r in use if r['seed'] in keep]
dens = sorted({r['dens'] for r in use})
# Keep only alphas whose whole density row has the full seed count (partially finished stages are left out).
cnt = {}
for r in use: cnt[(r['alpha'], r['dens'])] = cnt.get((r['alpha'], r['dens']), 0) + 1
target = max(cnt.values())
alphas = sorted(a for a in {k[0] for k in cnt} if all(cnt.get((a, d), 0) == target for d in dens))
use = [r for r in use if r['alpha'] in alphas]
lo = int(T * .75) + 100

def cut(r, key):
    ts = np.array([m['t'] for m in r['series'] if m['t'] <= T]); v = np.array([m[key] for m in r['series'] if m['t'] <= T])
    return ts, v

rows = []
for a in alphas:
    for d in dens:
        rs = [r for r in use if r['alpha'] == a and r['dens'] == d]
        if not rs: continue
        mean = np.mean([r['series'][0]['N'] for r in rs]) / A
        late = []
        for r in rs:
            ts, v = cut(r, 'clusterDensity'); late.append(v[ts >= lo].mean())
        _, frac = cut(rs[0], 'clusterFraction')
        lf = [cut(r, 'clusterFraction')[1][cut(r, 'clusterFraction')[0] >= lo].mean() for r in rs]
        rows.append({'alpha': a, 'dens_input': d, 'mean_density': round(mean, 4), 'seeds': len(rs), 'T': T,
                     'late_cluster_density': float(np.mean(late)), 'seed_min': float(min(late)), 'seed_max': float(max(late)),
                     'ratio_to_mean': float(np.mean(late) / mean), 'cutoff': math.ceil(2 * mean),
                     'late_cluster_mass_fraction': float(np.mean(lf))})
with open(out / f'summary-t{T}.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
(out / f'summary-t{T}.json').write_text(json.dumps(rows))

plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.labelcolor': '#19252d',
                     'text.color': '#19252d', 'xtick.color': '#19252d', 'ytick.color': '#19252d',
                     'figure.facecolor': '#fcfcfa', 'axes.facecolor': '#fcfcfa', 'savefig.facecolor': '#fcfcfa'})
cmap = plt.get_cmap('viridis')
col = {d: cmap(0.92 * i / max(1, len(dens) - 1)) for i, d in enumerate(dens)}
kfmt = FuncFormatter(lambda v, _: f'{v:g}' if v < 1000 else f'{v / 1000:g}k')
seedsTxt = sorted({r['seed'] for r in use}); nS = min(r['seeds'] for r in rows)

# 1. Small multiples: cluster density vs time, one panel per alpha, one curve per starting density (log-log).
def grid(alist, name):
  n = len(alist); ncol = 4; nrow = math.ceil(n / ncol)
  fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.55 * nrow + .6), layout='constrained', squeeze=False)
  for ax, a in zip(axes.flat, alist):
      for d in dens:
          rs = [r for r in use if r['alpha'] == a and r['dens'] == d]
          if not rs: continue
          ts, _ = cut(rs[0], 'clusterDensity'); v = np.mean([cut(r, 'clusterDensity')[1] for r in rs], axis=0)
          ax.plot(ts[1:], v[1:], color=col[d], lw=1.1)
      ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(25, T); ax.grid(alpha=.18, which='major')
      ax.xaxis.set_major_formatter(kfmt); ax.yaxis.set_major_formatter(kfmt); ax.yaxis.set_minor_formatter(NullFormatter())
      ax.set_title(f'α = {a:g}', loc='left', fontsize=10, fontweight='bold')
  for ax in axes.flat[n:]: ax.axis('off')
  for ax in axes[:, 0]: ax.set_ylabel('Cluster density\n(particles/site)')
  for ax in axes[-1, :]: ax.set_xlabel('Simulation step')
  sm = plt.cm.ScalarMappable(cmap=cmap, norm=LogNorm(6 * dens[0], 6 * dens[-1] / .92 * .92))
  cb = fig.colorbar(sm, ax=axes, shrink=.6, pad=.01); cb.set_label('Starting average density (particles/site)')
  fig.suptitle(f'Cluster density vs time · {len(dens)} starting densities · 90×90, centre + neighbours, sens 6 · '
               f'{"mean of " + str(nS) + " seeds" if nS > 1 else "seed " + str(seedsTxt[0])} · to {T:,} steps',
               x=.01, ha='left', fontsize=10, color='#5a676c')
  for ext in ['png', 'pdf']: fig.savefig(figdir / f'{name}-t{T}.{ext}', dpi=150)
  plt.close(fig)


sm = plt.cm.ScalarMappable(cmap=cmap, norm=LogNorm(6 * dens[0], 6 * dens[-1]))
base = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
grid([a for a in alphas if a in base], 'grid-time')
grid([a for a in alphas if a not in base or a in (.8, .9, 1)], 'grid-time-refine')
n = len(alphas)

# 2. Late level vs starting density, one line per alpha (log-log), with the 2x cutoff.
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), layout='constrained')
acm = plt.get_cmap('plasma')
for a in alphas:
    rr = [r for r in rows if r['alpha'] == a]; c = acm(0.9 * alphas.index(a) / max(1, n - 1))
    x = [r['mean_density'] for r in rr]
    axs[0].plot(x, [r['late_cluster_density'] for r in rr], '-o', ms=3, color=c, lw=1.3, label=f'α {a:g}')
    axs[1].plot(x, [r['ratio_to_mean'] for r in rr], '-o', ms=3, color=c, lw=1.3, label=f'α {a:g}')
xm = sorted({r['mean_density'] for r in rows})
axs[0].plot(xm, [math.ceil(2 * m) for m in xm], '--', color='#7a8488', lw=1, label='cutoff 2×avg')
axs[1].axhline(2, ls='--', color='#7a8488', lw=1)
for ax, yl in zip(axs, ['Late cluster density (particles/site)', 'Cluster density ÷ average density']):
    ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(alpha=.18); ax.set_xlabel('Average density (particles/site)'); ax.set_ylabel(yl)
    ax.xaxis.set_major_formatter(kfmt); ax.yaxis.set_major_formatter(kfmt)
axs[1].legend(frameon=False, fontsize=8, ncol=2, loc='upper right')
fig.suptitle(f'Late level (steps {lo:,}–{T:,}) · {nS} seed(s)', x=.01, ha='left', fontsize=10, color='#5a676c')
for ext in ['png', 'pdf']: fig.savefig(figdir / f'late-vs-density-t{T}.{ext}', dpi=150)
plt.close(fig)

# 3. Map: alpha x starting density, colour = cluster density / average density.
M = np.full((n, len(dens)), np.nan)
for r in rows: M[alphas.index(r['alpha']), dens.index(r['dens_input'])] = r['ratio_to_mean']
fig, ax = plt.subplots(figsize=(8.5, 0.34 * n + 1.8), layout='constrained')
im = ax.imshow(M, origin='lower', aspect='auto', cmap='magma', norm=LogNorm(2, np.nanmax(M)))
ax.set_xticks(range(len(dens)), [f'{6 * d:g}' for d in dens]); ax.set_yticks(range(n), [f'{a:g}' for a in alphas])
for i in range(n):
    for j in range(len(dens)):
        if not np.isnan(M[i, j]): ax.text(j, i, f'{M[i, j]:.0f}' if M[i, j] >= 10 else f'{M[i, j]:.1f}', ha='center', va='center', fontsize=6.5,
                                          color='white' if M[i, j] < 30 else 'black')
ax.set_xlabel('Starting average density (particles/site)'); ax.set_ylabel('α')
cb = fig.colorbar(im, ax=ax, pad=.01); cb.set_label('Cluster density ÷ average density')
ax.set_title(f'How many times denser than average the largest cluster is · steps {lo:,}–{T:,} · {nS} seed(s)', loc='left', fontsize=9.5, color='#5a676c')
for ext in ['png', 'pdf']: fig.savefig(figdir / f'map-ratio-t{T}.{ext}', dpi=150)
plt.close(fig)
# 4. Manik-style view: linear time, several alphas side by side, first 5k enlarged above the full run.
def focus(alist, name):
    alist = [a for a in alist if a in alphas]
    if not alist: return
    fig, axes = plt.subplots(2, len(alist), figsize=(4.1 * len(alist) + .8, 5.6), layout='constrained', squeeze=False,
                             gridspec_kw={'height_ratios': [2.2, 1]})
    for k, a in enumerate(alist):
        for d in dens:
            rs = [r for r in use if r['alpha'] == a and r['dens'] == d]
            if not rs: continue
            ts, _ = cut(rs[0], 'clusterDensity'); v = np.mean([cut(r, 'clusterDensity')[1] for r in rs], axis=0)
            for ax, stop in zip(axes[:, k], [5000, T]):
                u = ts <= stop; ax.plot(ts[u], v[u], color=col[d], lw=1.1)
        for ax, stop in zip(axes[:, k], [5000, T]):
            ax.set_yscale('log'); ax.set_xlim(0, stop); ax.grid(alpha=.18); ax.xaxis.set_major_formatter(kfmt)
            ax.yaxis.set_major_formatter(kfmt); ax.yaxis.set_minor_formatter(NullFormatter())
        axes[0, k].set_title(f'α = {a:g}', loc='left', fontsize=11, fontweight='bold')
        axes[1, k].set_xlabel('Simulation step')
    axes[0, 0].set_ylabel('Cluster density (particles/site)\nfirst 5,000 steps'); axes[1, 0].set_ylabel(f'full run')
    cb = fig.colorbar(sm, ax=axes, shrink=.7, pad=.01); cb.set_label('Starting average density (particles/site)')
    fig.suptitle(f'Cluster density vs time · 90×90, centre + neighbours, sens 6 · {"mean of " + str(nS) + " seeds" if nS > 1 else "seed " + str(seedsTxt[0])}',
                 x=.01, ha='left', fontsize=10, color='#5a676c')
    for ext in ['png', 'pdf']: fig.savefig(figdir / f'{name}-t{T}.{ext}', dpi=150)
    plt.close(fig)
focus([.8, .85, .9, 1], 'focus-transition'); focus([0, .4, .8], 'focus-low-alpha')

print(f'T={T}: {len(use)} runs, {n} alphas x {len(dens)} densities, seeds {seedsTxt}; wrote summary-t{T}.csv and figures.')
