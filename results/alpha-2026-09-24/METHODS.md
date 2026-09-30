# Alpha pilot — methods and reproducibility

Date: 2026-09-24. Purpose: prepare the next meeting with Manik, not estimate a final phase boundary. [interactive report](index.html).

## Rule and scope

Engine is extracted directly from `lgca-viz/lgca-lab.html`, `<script id="engineSrc">`. No independent dynamics implementation. Its SHA-256 is stored in `data/manifest.json`, every `run.json` and `validation.json`.

All pilot runs: boson occupancy, polar J=c_k·c_c, periodic axial hexagonal lattice, density parameter 0.4 per channel, WELL1024a seeds 12345 / 777 / 424242. Each run begins from its own seeded initialization. Runs of equal size/seed start from exactly the same lattice and generator state across alpha; diverging dynamics then consume the common stream according to their evolving occupancies. These are three replicate seeds, not 60 independent replicates of one condition.

The optional rule is field / M^alpha, with M the particle count in the same field stencil. Alpha=0 uses the original SUM arithmetic; alpha=1 uses the original AVERAGE arithmetic. Other alpha values divide max-subtracted boson scores before exponentiation. M=0 gives zero scores and uniform probabilities. The six-site neighbour ring excludes the centre. Own-site includes every particle at that site. There is no new radius control and no seven-site stencil.

The reference C++ is unchanged. The interpolation is our interpretation of the meeting notes; Zoom's `1/alpha` is recorded as an unresolved transcription issue, not implemented literally.

## Matrix: 60 runs, 348 saved states

| Group | Sizes L | Alpha | Seeds per condition | Sensitivity | Final step | Runs |
|---|---|---|---|---|---|---:|
| Main neighbours | 60, 120, 180 | 0, .5, 1 | 3 | 6 | 10,000 | 27 |
| Finer neighbours | 90 | 0, .25, .5, .75, .9, 1 | 3 | 6 | 10,000 | 18 |
| Own-site control | 60 | 0, .5, 1 | 3 | 6 | 10,000 | 9 |
| Matched mean-field control | 60 | 0, .5 | 3 | 6 / 14.4^(1−alpha) | 2000 | 6 |

The matched-field comparison reuses the main L=60, alpha=1 runs at t≤2000. Here 14.4=6 neighbours×6 channels×0.4 is the initial expected field population. This is an approximate scale control, not a pointwise normalization equivalence and not a way of holding all noise effects fixed.

Measurements every 50 steps. Exact states and RNG checkpoints at t=0,250,1000,2000,5000,10000 where within the run horizon. State is measured **after collision and streaming**, before the next collision. Large runs initially stopped at 2000, then resumed from their saved occupancy and all 33 WELL words to 10,000. An uninterrupted L=120 alpha=.5 seed=12345 run verifies the resumed final state and RNG exactly.

Initialization produces a fluctuating N around 2.4 L². N is measured and conserved within each run; ratios use each run's actual N. We do not replace it by 2.4 L².

## Stored data

- `data/<run-id>/state-t*.bin.gz`: gzip-compressed raw little-endian Int32 occupations. Offset `(j*L+i)*6+k`; no header.
- `data/<run-id>/run.json`: settings, engine hash, complete sampled time series, saved-state hashes, RNG states, occupation histograms, threshold checks and boundary profiles.
- `data/manifest.json`: planned settings/horizons. The builder refuses incomplete runs.
- `measurements.csv`: time series with rule, size, seed and sensitivity attached to every row.
- `summary.json`: joined metadata and measurements, also embedded in the report so no data server/API is needed.
- `tool-settings.json`: 60 recipes importable into Lab. Applying one resets the simulation; it does not load a late-time state.
- `validation.json`: data validation and one-box-crossing recurrence checks.

## Definitions

For n_x=sum_k n_xk and N=sum_x n_x:

- **nmax** = max_x n_x; **kmax** = max_{x,k} n_xk. Site versus site-channel peaks are intentionally kept distinct.
- **siteFraction** = nmax/N; **channelFraction** = kmax/N.
- **effectiveArea** = N²/sum_x n_x². This participation area refers to all particles, and is not a geometric cluster area.
- **polar** = norm(sum_{x,k} n_xk c_k)/N.
- **localPolar** = average over occupied sites of norm(sum_k n_xk c_k)/n_x. A singleton has localPolar=1 regardless of global order.
- **weightedLocalPolar** = sum_x norm(sum_k n_xk c_k)/N; particle weighting differs from site weighting.
- **histogram** = number of sites with each occupation n, including zeros.

**Clusters:** retain sites n_x≥ceil(f N/L²), using f=2 for the main analysis. Take connected components of the six-neighbour graph with periodic wrapping. Select the component with greatest particle mass; exact ties follow site traversal order. Report its mass, number of sites (area), mass/area, mass/N, and vector polarization. Threshold checks f=1.5 and 4 are saved at each exact checkpoint. We do not smooth before thresholding. With the mean near 2.4, the main threshold is 5, and the auxiliary thresholds are 4 and 10.

The largest component need not be the one containing nmax. Identity is selected independently at each measurement; a jump in cluster metrics may mean a different component has become largest. At t=0, connected dense fluctuations should not be interpreted as established physical clusters. A threshold-defined component may wrap around the periodic box; no Euclidean radius or planar convex hull is inferred.

**Edge profiles:** select the largest f=2 component. Boundary sites are its sites adjacent to any site outside that component. Multi-source breadth-first search on the full periodic graph measures the minimum hop distance from those boundary sites. Boundary =0, interior distances negative, exterior positive. Report bins up to 12 hops, including bin site/particle counts, mean occupation across all sites in the bin, mean local order across occupied sites, and particle-weighted order. No arbitrary order is assigned to empty bins. Other clusters may occur in exterior bins. The illustrated 180×180 profiles average equally over the three seeds only at distances where all three contribute; for order each must have occupied sites in the bin. Unavailable depths are absent. Zero density is omitted from the logarithmic plot.

This is a distance-to-boundary profile, not a radial profile from a centroid, not a moving-cluster time average, and not a front-versus-rear decomposition. Those alternatives require an agreed definition and additional tracking.

## Figures and uncertainty

Density uses the Lab's palette with fixed index `floor(255 log(1+n)/log(1001))`, clipped to [0,255], empty sites dark. Same scale for every panel; no per-frame autoscaling. Hexagonal shear: site (i,j) maps to x=(2i+j) mod 2L, y=2j, with a 2×2 block. Display height is scaled by sqrt(3)/2. All figures use exact saved states. Direction uses the existing tool renderer. Local order uses a linear 0–1 palette.

Ensemble tables/figures use arithmetic means and **sample standard deviations (ddof=1)** across three seeds. Error bars are not confidence intervals. The displayed mean of each run's mass/area is not generally the ratio of mean mass to mean area. The static figures are independent of the report's selectors and are explicitly labeled.

## Validation completed

- Required reference suite: **All 10 C++/JavaScript parity scenarios passed.** See `parity.log`.
- Alpha endpoints: 16 legacy/endpoint trajectories, 120 steps each, exact full-state and all-word RNG equality across both models and all existing boson toggles.
- Intermediate alpha: conservation, exclusion for fermions, zero-sensitivity equality, independent asymmetric-field CDF, isolated empty-neighbour field.
- Real C++ adversarial checks: 24 state/RNG checkpoints across four 300-step runs, 96×60, sens=8, boson dens=.4/2.25, fermion dens=.2/.95. See `alpha-validation.log`.
- Data: all 60 runs and 348 checkpoints have matching engine/state hashes, nonnegative counts, conserved N, and independently recomputed recorded statistics. Periodic cluster fixture verifies mass, area, boundary density and wrap connectivity. Uninterrupted/resumed state+RNG equality verified.
- At t=10,000, all 45 neighbour runs were advanced one lattice-side duration L. **6/45** reproduce their entire occupancy array exactly; see the individual IDs in `validation.json`. The RNG continues advancing, so these are observed state recurrences, not a proof of permanent periodicity of the full stochastic state. Three recurrences occur for alpha=0 and three for alpha=.5 at L=60/90.

## Limits of interpretation

No stationary distributions, universal alpha_c, asymptotic exponents, or infinite-time divergence have been established. Finite-time plateaus can reflect long transients, weak mixing of ballistic structures on a periodic lattice, extremely small escape probabilities or floating-point saturation. The C++ and JS share the finite-precision limitation: sufficiently separated softmax weights can produce CDF entries rounded to one. Exact parity establishes reference fidelity, not rare-event accuracy at arbitrarily long times.

The important next choices are the field stencil, definition of condensation, cluster definition and relaxation protocol. Our size scan at one fixed time is a diagnostic, not sufficient finite-size scaling evidence. Our alpha scan explores one sensitivity/density slice; a crossover bracket does not define a universal phase boundary.

## Reproduce

From the repository root, with Node and C++ on PATH:

```sh
node lgca-viz/verify-parity.mjs
node lgca-viz/verify-alpha.mjs
node results/alpha-2026-09-24/build/run.mjs
node results/alpha-2026-09-24/build/validate.mjs
node results/alpha-2026-09-24/build/package.mjs
python3 results/alpha-2026-09-24/build/plots.py
```

The plot script requires numpy and matplotlib. Runs use only Node's standard library and the shipped engine. Completed jobs are reused; shorter completed jobs resume from checkpoints. The runner rejects a changed engine when resuming. To deliberately rerun from scratch, use a new output directory or explicitly move the old data out first. It runs up to three independent OS processes; the general simulator's architecture is unchanged.

The alpha test uses the pre-change engine snapshot in `build/lgca-lab.before.html` as its endpoint regression baseline. This is a historical fixture, not a second active version of the tool.

Open `index.html` directly, or serve the repository locally; keep its `images/` and `figures/` folders next to it. The link to the general tool assumes the current repository layout.
