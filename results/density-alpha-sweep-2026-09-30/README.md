# Density × alpha sweep (Sep 30)

Polar bosons, field from centre site + 6 neighbours, sensitivity 6, field divided by
M^alpha.

- Main sweep: box 90×90, 12 starting densities (6 × dens = 0.15 to 14.4 particles per
  site), 25 alpha values from 0 to 1, 5 seeds (3 for the finer alphas), 20,000 steps.
- Box size runs: dens 0.4, L = 90, 120, 180, 240, seeds 12345, 777 and 424242.

| File | What is inside |
|---|---|
| `figures/focus-transition-t20000-prelim.png` | Cluster density vs time for alpha 0.8, 0.85, 0.9, 1. Top row: time linear until step 500, log after. Bottom row: all the run, linear |
| `figures/grid-time-t20000-prelim.png` | Same curves for alpha 0 to 1 every 0.1 (log–log) |
| `figures/grid-time-refine-t20000-prelim.png` | Same curves for alpha 0.75 to 1, more fine |
| `figures/condensate-vs-box-size-t20000.png` | Late condensate size N_c vs box side |
| `cluster-density-vs-time-t20000.csv` | Cluster density (mean over seeds) at each saved step, one column for each alpha and dens, e.g. `a0.85_d0.4` |
| `summary-t20000-prelim.csv` | Late averages (steps 15,100 to 20,000) for each alpha and dens, with number of seeds |
| `condensate-vs-box-size-t20000.csv` | Late N_c for each alpha and box size, mean and min/max over seeds |

`python3 build/figures.py` makes the figures from the CSVs (needs NumPy and
Matplotlib). With `--data <folder>` it first makes the CSVs again from the saved runs,
these are not in the repo because they are too big.
