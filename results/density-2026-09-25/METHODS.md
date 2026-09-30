# Density sweep — 25 September 2026

Polar bosons; centre plus six neighbours (`both`); power normalization by the total interacting particle count to exponent alpha. Sensitivity 6, periodic 90×90 box. Density inputs 0.1, 0.2, 0.4, 0.8 per channel (expected site means 0.6, 1.2, 2.4, 4.8); alpha 0.8, 0.9, 1; seeds 12345, 777, 424242, 20260925, 314159. These 60 fresh trajectories run through 80,000 steps each. No particle injection. Actual particle totals are measured and conserved within each trajectory.

This varies mean density at fixed box size. The previous neighbours-only size sweep varies box size at fixed mean density. Both interaction field and comparison axis differ: their differences cannot be attributed to either change alone. A matched interaction comparison at density input0.4 and size90 is possible for alpha0.8/0.9 and the shared five seeds.

Cluster = largest-mass periodic six-connected component of sites with at least ceil(2 × actual mean occupancy). Density = cluster mass / area in lattice sites. Alternative thresholds1.5×/4× are retained at checkpoints; maximum site/channel occupancy fractions and mass/area are retained separately. The selected largest cluster can change identity. Zero cluster area is recorded as zero density, not missing data. Because the threshold scales with system mean, interpreting a density sweep requires those threshold checks.

Observables are sampled every25 steps through5000 and every100 thereafter. Exact Int32 lattice states and all33 RNG words are saved at0,100,500,1000,2000,5000,10000,20000,40000,60000,80000. Each checkpoint has SHA-256 and count-conservation checks. The display emphasizes0–5000 steps and keeps the full horizon visible underneath with the same density scale. This is only a visual zoom; no late evolution is omitted from computation or analysis. Curves show instantaneous means ± sample SD over five seeds. Late numbers average each trajectory over the final quarter, then summarize five seed means. No asymptote is fitted or asserted from a visually flat segment.

The pinned engine source and hash are in build/engine-source.js and provenance.json. All core tests passed before the sweep; the centre-inclusive variant has independent analytical CDF checks. Completed-run validation remeasures every stored state, verifies hashes/conservation/RNG shape, paired initializations across alpha, and two independent exact state+RNG continuations.

From repository root:

```sh
node results/density-2026-09-25/build/run.mjs
node results/density-2026-09-25/build/analyze.mjs
node results/density-2026-09-25/build/validate.mjs
python3 results/density-2026-09-25/build/plots.py
node results/density-2026-09-25/build/package.mjs
node results/alpha-refinement-2026-09-25/build/package.mjs
```

The same results page displays both experiments; no separate competing report. All PNG/PDF plots are under figures/. The runner checkpoints atomically and resumes exact states, with three independent worker processes and no engine optimization.
