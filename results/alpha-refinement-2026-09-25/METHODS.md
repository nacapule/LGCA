# Alpha refinement — 25 September 2026

Question: on either side of the candidate crossover, does attained cluster density increase with total particles at fixed mean density, or approach a size-independent value?

## Controlled comparison

- Bosons, polar alignment, six neighbours excluding the centre, periodic boundaries.
- Divide the polar field by M^α. Sensitivity 6, density input 0.4 per channel (expected mean 2.4 particles/site).
- α: 0.75, 0.80, 0.85, 0.90. Sizes: 90×90, 120×120, 180×180.
- Five seeds per condition: 12345, 777, 424242, 20260925, 314159. 60 runs total.
- Fresh initialization for every run; matching size/seed pairs have exactly the same initial lattice and RNG across α. No particle injection.
- Initial horizon: 20,000 steps. Metrics every 100; exact occupancy and all 33 WELL RNG words at 0, 1000, 5000, 10000 and 20000. If the time comparison shows continuing change, extend the sweep to 40000, then at most 80000 in this exploratory study. Store final horizons explicitly in the manifest.

Density stays fixed while size changes. Doubling L approximately quadruples N. Changing the density input is a separate parameter study, not a substitute for this comparison.

## Primary measurement

Cluster density = cluster mass / cluster area, in particles per lattice site. Retain sites with count at least ceil(2 × the run's measured mean). Take connected components on the periodic six-neighbour graph and select the component with largest mass. Area is its number of sites. Also record thresholds at 1.5× and 4× the mean, because cluster boundaries depend on that choice.

Record mass and area separately, largest-site and largest-channel counts divided by N, and participation area. A large mass fraction alone does not prove growing internal density. The selected largest component can change identity; its time series does not follow one named cluster.

The main table and size plots average each run over the last quarter up to the selected horizon, then show the five-seed mean and sample standard deviation of those per-run averages. At horizon 40000, that uses samples 30100…40000. Time-curve plots show instantaneous five-seed means and standard deviations. Snapshot metrics remain instantaneous; their numbers need not match the time-averaged table. The average of mass/area is not generally mean mass divided by mean area. Consecutive observations are correlated; they are not independent replicates. Compare time windows in the penultimate and final quarters of each run. Relative drift is a descriptive warning of continuing evolution, not a stationarity test. Small drift cannot rule out slow coarsening or metastability.

Evidence sought: as L increases, cluster density rises while a nonvanishing fraction of N remains concentrated in a subextensive area. A roughly constant density with area growing with mass instead suggests an extended finite-density phase. Three sizes and a finite observation horizon do not establish an infinite-system limit or a unique critical alpha. There can also be multiple long-lived clusters whose masses fail to scale with N over these sizes and times.

## Injection thought experiment

Continuous injection into a fixed box changes N/area and the dynamical protocol. Even a uniform state must grow denser if injection continues without removal. Stop-and-relax injection could explore history dependence or serve as an alternative initialization, but would need comparison with fresh runs at the same final settings and sufficient relaxation. For an injection study, normalized mass fractions or density relative to the rising system mean would be more informative than absolute density alone. That probes a different limit and cannot by itself replace increasing both N and system area at fixed mean density.

The original own-site condensate study defines a macroscopic occupation in a large-system limit; our neighbour rule and extended threshold-defined clusters require their own checks: [Nava-Sedeño et al., polar condensates and nematic filaments](https://arxiv.org/html/2402.04450v1).

## Reproducibility

The study originally extracted the engine directly from `lgca-viz/lgca-lab.html`; every run stores its SHA-256. That exact source is now pinned at `../alpha-2026-09-24/build/engine-source.js`, and the runner/validator load it. Its hash remains identical to all recorded runs. This preserves reproduction after the Lab gained the separate centre-plus-neighbours option. Analysis reuses the validated periodic-cluster measurement code from the September 24 pilot. Counts and RNG checkpoint exactly after streaming. Checkpoint JSON is written atomically; resuming validates settings and state hash.

From the repository root:

```sh
node results/alpha-refinement-2026-09-25/build/run.mjs
node results/alpha-refinement-2026-09-25/build/analyze.mjs
node results/alpha-refinement-2026-09-25/build/validate.mjs
python3 results/alpha-refinement-2026-09-25/build/plots.py
node results/alpha-refinement-2026-09-25/build/package.mjs
```

The plot script needs numpy and matplotlib. Findings are an explicit interpretation in `findings.json`, tied to the SHA-256 of `analysis.json`; after changing the data, reassess them and update that hash before rebuilding the page.

The completed study ends at 80000, now the runner default. The recorded sequence was 20000 → 40000 → 80000, always resuming exact states. The runner refuses to lower the recorded horizon. Use a separate output directory for an independent shorter study. The small worker pool runs separate experiments; no simulator architecture or engine optimization is involved.

## Why an apparent intermediate threshold may move

An elementary retention calculation also applies to a specially prepared neighbour configuration. Put n equally oriented particles on each of two adjacent sites, with all other sites empty. Each site sees the other site's n particles; moving together preserves their adjacency. The forward-channel probability is

`p = 1 / [1 + 2 exp(−u/2) + 2 exp(−3u/2) + exp(−2u)]`, with `u = 6 n^(1−α)`.

The expected number changing direction across both sites is `2n(1−p)`. For every fixed α < 1 it tends to zero as n grows; at α = 1 the departing fraction stays positive. This is a calculation for a prepared, aligned pair, not proof that a random initial state forms a thermodynamic condensate. For large u, the expected number changing direction is approximately `4n exp[−3 n^(1−α)]`. Its turning point is `n* ≈ [1 / (3(1−α))]^(1/(1−α))`: about 13 particles per site at α=.80, 205 at .85, and 169,000 at .90. These are turning points of this particular escape proxy, not cluster-formation thresholds or measured critical densities. They illustrate how sharply the accessible mass scale can change with alpha, and motivate asking Manik whether the apparent intermediate boundary could drift toward the endpoint α=1. It explains why a sharp finite-size crossover strictly below 1 should not automatically be declared a critical alpha. Formation, mergers, fragmentation and stationary mass distributions must also be understood.

## Files and validation

`data/<run-id>/run.json` holds parameters, engine hash, sampled measurements, checkpoint hashes and RNG words. `state-t*.bin.gz` is gzip-compressed little-endian Int32 occupancy, indexed as `(j*L+i)*6+k`. These are after-streaming states. The initial N fluctuates around 2.4 L² according to the reference initialization; measurements use actual N.

`analysis.json` contains per-condition means, seed spreads, late time windows, descriptive drift and threshold checks. `measurements.csv` retains every sampled trajectory. Earlier 20,000- and 40,000-step analyses are preserved separately. The browser page embeds its numerical summaries and uses local PNGs; it also works offline when the folder structure is retained.

Validation checks every checkpoint hash, nonnegative counts, conservation, recomputed cluster measurements at all three thresholds, identical initialization across alpha at equal size/seed, and all six shared pilot trajectories at step 10000. A separately reconstructed 120×120, alpha .85, seed12345 continuation from step10000 must match the saved final state and all RNG words exactly. The engine hash must match the engine already verified against the reference C++ in the September24 pilot.

A temporal plateau does not uniquely identify its cause. Separated coherently moving clusters may fail to meet over the observation window; escape probabilities can be very small; finite-precision CDFs and the 32-bit uniform RNG limit extremely rare transitions. Exact agreement with the reference engine establishes reproducibility, not accurate sampling of arbitrarily rare events or equilibration.
