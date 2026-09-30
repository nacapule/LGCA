# Density × alpha sweep (September 30)

Polar bosons, field from the centre site plus its six neighbours, sensitivity 6,
field divided by M^alpha. Main sweep: 90×90 box, 12 starting densities
(6 × dens = 0.15 to 14.4 particles/site), 25 values of alpha from 0 to 1,
5 seeds (3 for the finer alphas), 20,000 steps. Box-size runs: dens 0.4,
L = 90, 120, 180, 240, seeds 12345, 777 and 424242.

| File | Contents |
|---|---|
| `figures/focus-transition-t20000-prelim.png` | Cluster density vs time at alpha 0.8, 0.85, 0.9, 1. Top: time axis linear up to step 500 and logarithmic after, so the early rise is wide; bottom: full run on a linear axis |
| `figures/grid-time-t20000-prelim.png` | The same curves for alpha 0 to 1 in steps of 0.1 (log–log axes) |
| `figures/grid-time-refine-t20000-prelim.png` | The same curves for alpha 0.75 to 1 in finer steps |
| `figures/condensate-vs-box-size-t20000.png` | Late condensate size N_c against box side |
| `cluster-density-vs-time-t20000.csv` | Seed-mean cluster density at each saved step; one column per alpha and dens (for example `a0.85_d0.4`) |
| `summary-t20000-prelim.csv` | Late-window averages (steps 15,100–20,000) for every alpha and dens, with seed counts |
| `condensate-vs-box-size-t20000.csv` | Late N_c for each alpha and box size: seed mean and range |

`python3 build/figures.py` redraws the figures from the CSV files (needs NumPy and
Matplotlib). With `--data <folder>` it first recomputes the CSV files from the
saved run records, which are not included here because of their size.
