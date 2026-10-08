# LGCA

A hexagonal lattice-gas cellular automaton with a C++ simulation and a browser
Lab for studying alignment and clustering.

## Open the Lab and the island study

- **[Open LGCA Lab](https://nacapule.github.io/LGCA/lgca-viz/lgca-lab.html)**
- **[Open the island study](https://nacapule.github.io/LGCA/results/islands-2026-10-06/island-study.html)**
  (Oct 7, 2026) — island density over time for each particle count, with two
  ways of adding particles and three ways of choosing the island.

The Lab supports fermion and boson models, adjustable density and alignment,
replay, site inspection, density views and measurement exports. Settings
can be saved and shared.

## Island study

This study changes two things in the density and alpha study below. The
island is a whole group of connected occupied sites, with no density
cutoff. And
the particle count N grows in two ways: a denser 90 × 90 box, as before,
or a larger box at a fixed density. The page lays the old and new versions
of the same plot over each other.

| Data | Contents |
|---|---|
| [Island density over time](results/islands-2026-10-06/island-density-vs-time.csv) | Seed means for each box, alpha and N: `old` is the cluster of the earlier study, `new` the island of the busiest site, `densest` the island with the most particles per site |
| [Late island density](results/islands-2026-10-06/late-island-density.csv) | Late means for three choices of island, with seed ranges |
| [Islands of the saved lattices](results/islands-2026-10-06/islands-at-saved-steps.csv) | Particles and sites of the busiest site's island, the heaviest island and the densest island in each saved lattice |

## Density and alpha study

These runs compare cluster density at different starting densities and
values of **alpha**. Alpha 0 uses the alignment field as a sum; alpha 1
uses its average per particle. The runs use polar bosons, a field from the
central site and six neighbours, and sensitivity 6.

| Figure | Contents |
|---|---|
| **[Cluster density over time](results/density-alpha-sweep-2026-09-30/figures/focus-transition-t20000-prelim.png)** | Alpha 0.8, 0.85, 0.9 and 1 at 12 starting densities; early growth above, full runs below |
| [Coarse alpha comparison](results/density-alpha-sweep-2026-09-30/figures/grid-time-t20000-prelim.png) | Alpha 0 to 1 in steps of 0.1 |
| [Finer alpha comparison](results/density-alpha-sweep-2026-09-30/figures/grid-time-refine-t20000-prelim.png) | More closely spaced alphas between 0.75 and 1 |
| [Condensate occupation across box sizes](results/density-alpha-sweep-2026-09-30/figures/condensate-vs-box-size-t20000.png) | Maximum occupation of one channel at one site, N_c, for box sides 90, 120, 180 and 240 |

The figures cover runs through 20,000 steps, averaged across several seeds.

| Data | Contents |
|---|---|
| [Cluster-density time curves](results/density-alpha-sweep-2026-09-30/cluster-density-vs-time-t20000.csv) | Mean curves for each alpha and starting density |
| [Late cluster-density averages](results/density-alpha-sweep-2026-09-30/summary-t20000-prelim.csv) | Averages over steps 15,100–20,000, with seed counts and ranges |
| [Condensate occupation by box size](results/density-alpha-sweep-2026-09-30/condensate-vs-box-size-t20000.csv) | Late N_c means and ranges across seeds |

## Documentation and source

| Document | Contents |
|---|---|
| [Model and observables](docs/PHYSICS.md) | Collision rules, field settings, alpha, measurement and island definitions, island CSV columns |
| [Lab C++ source](native/README.md) | Simulation source and instructions for rebuilding after code changes |
| [Reference C++](lgca/lgca_clean-1.cpp) | The reference model implementation |
| [Earlier C++ variant](lgca/lgca_noib.cpp) | Neighbour fields, polar bosons and normalization |
| [Original repository](https://github.com/tailswalker/LGCA) | tailswalker/LGCA, which this repository is forked from |

## Resumen

LGCA Lab permite explorar el modelo en el navegador, sin instalar nada:
cambiar parámetros, pausar, avanzar o retroceder, mirar la densidad y exportar
mediciones. La simulación corre en C++.

El estudio compara la densidad de los cúmulos para distintas densidades
iniciales y valores de alpha, y la ocupación máxima de un canal al aumentar
el tamaño de la caja. El estudio de islas (7 de octubre de 2026) cambia dos
cosas: la isla es un conjunto conectado de sitios ocupados, sin umbral de
densidad, y el número de partículas también crece con una caja más grande a
densidad fija, no solo con una caja de 90 × 90 más densa. Su página permite comparar las versiones de cada gráfica.
Las tablas de arriba enlazan las figuras, los datos y las definiciones del
modelo.
