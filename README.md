# LGCA

This is a fork of the LGCA code (hexagonal lattice gas, with fermion and boson versions).
The original C++ is `lgca/lgca_clean-1.cpp`, with the WELL1024a random generator in
`rng/` and CImg in `CImg/` for the window display. We did not change these files.

What we add is two things: a version of the model that runs in the browser (LGCA Lab),
and the results of the alpha study we talked about in September.

## LGCA Lab

Open it here: https://nacapule.github.io/LGCA/lgca-viz/lgca-lab.html, or open
`lgca-viz/lgca-lab.html` directly.

The Lab is the C++ translated to JavaScript, random generator included, so with the
same seed and the same parameters it gives the same lattice in every step and the same
random numbers than the C++. Default parameters are the ones of the C++. The alpha
normalization and the field with centre site + 6 neighbours are options added for
the study.

To check it against the C++: `node lgca-viz/verify-parity.mjs` (needs Node.js and a
C++ compiler). It compiles the C++ and compares the whole lattice and the full state of
the random generator with the JavaScript, in 10 cases with fermions and bosons,
different seeds, densities and sensitivities. More detail in [docs/FIDELITY.md](docs/FIDELITY.md).

## Alpha study (Sep 30)

Question: for which alpha the cluster density keeps growing when we put more particles.
The alignment field is divided by M^alpha, with M the number of particles that
contribute to the field. Alpha 0 is the original sum, alpha 1 is the average. We use
polar bosons, field from the centre site and its 6 neighbours, sensitivity 6.
Definitions in [docs/PHYSICS.md](docs/PHYSICS.md).

- [Cluster density vs time for alpha 0.8, 0.85, 0.9, 1](results/density-alpha-sweep-2026-09-30/figures/focus-transition-t20000-prelim.png).
  Each curve is one starting density (12 of them, 0.15 to 14.4 particles per site),
  box 90×90. In the top row the time axis is linear until step 500 and log after, to
  see the first growth.
- Same curves for all the alphas: [0 to 1](results/density-alpha-sweep-2026-09-30/figures/grid-time-t20000-prelim.png),
  and [0.75 to 1 more fine](results/density-alpha-sweep-2026-09-30/figures/grid-time-refine-t20000-prelim.png).
- [Condensate size vs box size](results/density-alpha-sweep-2026-09-30/figures/condensate-vs-box-size-t20000.png),
  density fixed at 2.4 particles per site, box side 90, 120, 180, 240. N_c is the
  max number of particles in one channel of one site, like in the 2023 paper. The late
  value goes 91 → 170 → 700 → 1887 for alpha 0.9, and 24 → 30 → 40 → 69 for alpha 1.

Curves are the mean over seeds (5, or 3 for the finer alphas), runs of 20,000 steps.
The numbers of each figure are in CSV files in the [study folder](results/density-alpha-sweep-2026-09-30/),
and `build/figures.py` makes the figures again from them.

## Resumen

LGCA Lab es el código en C++ pasado a JavaScript; con la misma semilla da la misma red
en cada paso y los mismos números aleatorios (lo comprobamos en 10 casos). Se abre en el
navegador, sin instalar nada. El estudio más reciente busca para qué alpha la
densidad del cúmulo sigue creciendo al meter más partículas: las figuras son la
densidad del cúmulo contra el tiempo para 12 densidades iniciales, y el tamaño del
condensado N_c contra el tamaño de la caja.
