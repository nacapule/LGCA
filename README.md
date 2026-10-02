# LGCA

A hexagonal lattice-gas cellular automaton with a C++ simulation and a browser
Lab for studying alignment and clustering.

## Open the Lab

- **[Open LGCA Lab](https://nacapule.github.io/LGCA/lgca-viz/lgca-lab.html)**
- **[Lab HTML file](lgca-viz/lgca-lab.html)** — download this file and open it
  in a browser to work offline. Nothing to install.

The Lab supports fermion and boson models, adjustable density and alignment,
replay, site inspection, density views and measurement exports. Settings
can be saved and shared.

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
Full settings and measurement definitions are in the study guide below.

| Data | Contents |
|---|---|
| [Cluster-density time curves](results/density-alpha-sweep-2026-09-30/cluster-density-vs-time-t20000.csv) | Mean curves for each alpha and starting density |
| [Late cluster-density averages](results/density-alpha-sweep-2026-09-30/summary-t20000-prelim.csv) | Averages over steps 15,100–20,000, with seed counts and ranges |
| [Condensate occupation by box size](results/density-alpha-sweep-2026-09-30/condensate-vs-box-size-t20000.csv) | Late N_c means and ranges across seeds |

## Documentation and source

| Document | Contents |
|---|---|
| [Model and observables](docs/PHYSICS.md) | Collision rules, field settings, alpha and measurement definitions |
| [Study guide](results/density-alpha-sweep-2026-09-30/README.md) | Run settings, CSV columns and how to redraw the figures from the data |
| [Lab C++ source](native/README.md) | Simulation source and instructions for rebuilding after code changes |
| [Reference C++](lgca/lgca_clean-1.cpp) | The original model implementation |
| [Earlier C++ variant](lgca/lgca_noib.cpp) | Neighbour fields, polar bosons and normalization |
| [Original repository](https://github.com/tailswalker/LGCA) | tailswalker/LGCA, which this repository is forked from |

## Resumen

LGCA Lab permite explorar el modelo en el navegador, sin instalar nada:
cambiar parámetros, pausar, avanzar o retroceder, mirar la densidad y exportar
mediciones. La simulación corre en C++.

El estudio compara la densidad de los cúmulos para distintas densidades
iniciales y valores de alpha, y la ocupación máxima de un canal al aumentar
el tamaño de la caja. Las tablas de arriba enlazan las figuras, los datos
y las definiciones del modelo.
