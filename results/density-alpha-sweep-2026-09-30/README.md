# Density and alpha study

These runs ask how cluster density changes with the number of available
particles and with alpha. They use polar bosons, a field from the centre
and six neighbours, sensitivity 6, and division by M^alpha.

## Figures

| Figure | Contents |
|---|---|
| **[Cluster density over time](figures/focus-transition-t20000-prelim.png)** | Alpha 0.8, 0.85, 0.9 and 1, with a curve for each starting density. Upper row: time is linear to step 500 and logarithmic after it. Lower row: the full run on a linear time axis. |
| [Coarse alpha comparison](figures/grid-time-t20000-prelim.png) | Alpha 0 to 1 in steps of 0.1, on log–log axes |
| [Finer alpha comparison](figures/grid-time-refine-t20000-prelim.png) | More closely spaced alphas from 0.75 to 1 |
| [Condensate occupation across box sizes](figures/condensate-vs-box-size-t20000.png) | Late N_c across box sides 90, 120, 180 and 240, with the range across seeds |

## Settings

| Comparison | Settings |
|---|---|
| Starting density | 90 × 90 box; 12 expected densities from 0.15 to 14.4 particles/site; 25 alphas from 0 to 1 |
| Box size | Starting density input 0.4, about 2.4 particles/site; sides 90, 120, 180 and 240; alphas 0.8, 0.85, 0.9 and 1 |
| Seeds | Five for the coarse density curves, three for the finer curves and box-size comparison |
| Duration | 20,000 steps per included run |
| Late averages | Steps 15,100–20,000 within each run, then averaged across seeds |

Cluster density is the mass divided by the area of the largest-mass
connected component above the density threshold. N_c is the maximum
occupation of one channel at one site. These measure different things.

The runs show finite-size and finite-time behavior. Their late levels
are measured averages, not fitted infinite-time limits.
[Model and measurement definitions](../../docs/PHYSICS.md).

## CSV data

| File | Contents |
|---|---|
| [Time curves](cluster-density-vs-time-t20000.csv) | Mean cluster density at each sampled step, one column per alpha and density input |
| [Late cluster-density averages](summary-t20000-prelim.csv) | Means, number of seeds, and smallest/largest seed means |
| [Condensate occupation by box size](condensate-vs-box-size-t20000.csv) | Late N_c means and seed ranges |

A time-curve column such as `a0.85_d0.4` means alpha 0.85 and density
input 0.4 per channel. The expected starting density is six times that input.

## Redraw the figures from the CSV data

This command redraws the four figures from the included CSVs. It does not
run new simulations. It needs Python, NumPy and Matplotlib.

From the repository folder:

```sh
python3 results/density-alpha-sweep-2026-09-30/build/figures.py
```

This reads the included CSVs and writes the four PNGs. To rebuild the CSVs
from saved runs before plotting:

```sh
python3 results/density-alpha-sweep-2026-09-30/build/figures.py --data /path/to/saved/runs
```

The raw checkpoint states are kept in the working research checkout and
are not included in this public copy.
