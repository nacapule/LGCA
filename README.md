# LGCA

A hexagonal lattice-gas cellular automaton for studying how particles align
and cluster. The C++ simulation runs directly in the browser.

- **[Open LGCA Lab](https://nacapule.github.io/LGCA/lgca-viz/lgca-lab.html)** —
  change the rules, replay a run, inspect sites and export measurements.
- **[Open the island study dashboard](https://nacapule.github.io/LGCA/results/islands-2026-11-17/island-study.html)** —
  explore island density and spatial entropy across box sizes, α and time.

## Island study

An island is a connected group of occupied sites. The dashboard compares three ways
to choose an island. The runs add particles either by filling a fixed box more
densely or by enlarging the box at a fixed density.

The larger-box runs reach 1.28 million steps for every α tested except α = 0, which
reaches 320,000. Use **vs L** to compare boxes at the same step and **single seeds**
to see how individual runs differ.

[Island density over time](results/islands-2026-11-17/island-density-vs-time.csv) ·
[Late densities and seed ranges](results/islands-2026-11-17/late-island-density.csv) ·
[Measurements from saved lattices](results/islands-2026-11-17/islands-at-saved-steps.csv) ·
[Spatial entropy data](results/islands-2026-11-17/entropy-page.json)

The [October island study](https://nacapule.github.io/LGCA/results/islands-2026-10-06/island-study.html)
is available with [its data](results/islands-2026-10-06/).

## Density and α study

These runs use polar bosons with alignment sensitivity β = 6. The field includes
each site and its six neighbours. At α = 0 the field sums their contributions.
At α = 1 it averages them per particle.

- **[Cluster density over time](results/density-alpha-sweep-2026-09-30/figures/focus-transition-t20000-prelim.png)** —
  α = 0.8, 0.85, 0.9 and 1 across 12 starting densities, through 20,000 steps.
- **[Coarse](results/density-alpha-sweep-2026-09-30/figures/grid-time-t20000-prelim.png)
  and [fine](results/density-alpha-sweep-2026-09-30/figures/grid-time-refine-t20000-prelim.png) α comparisons** —
  α = 0 to 1 in steps of 0.1, and a closer look from 0.75 to 1.
- **[Condensate occupation by box size](results/density-alpha-sweep-2026-09-30/figures/condensate-vs-box-size-t20000.png)** —
  maximum occupation of a single channel at one site, for box sides 90 to 240.

[Cluster density time series](results/density-alpha-sweep-2026-09-30/cluster-density-vs-time-t20000.csv) ·
[Late cluster densities](results/density-alpha-sweep-2026-09-30/summary-t20000-prelim.csv) ·
[Condensate occupation data](results/density-alpha-sweep-2026-09-30/condensate-vs-box-size-t20000.csv)

## Model and source

- [Model and measurements](docs/PHYSICS.md) covers the collision rules, alignment and observables.
- [C++ simulation and command-line tools](native/README.md) covers research runs and rebuilding the Lab.
- [Reference C++](lgca/lgca_clean-1.cpp) · [Earlier C++ variant](lgca/lgca_noib.cpp) ·
  [Original repository](https://github.com/tailswalker/LGCA)

## Resumen

El Lab permite cambiar las reglas, inspeccionar sitios y exportar mediciones desde
el navegador. El estudio de islas muestra la densidad de las islas y la entropía
espacial a lo largo del tiempo y para distintos tamaños de caja.
