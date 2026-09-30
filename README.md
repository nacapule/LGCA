# LGCA

Hexagonal lattice-gas cellular automaton. The reference code is
`lgca/lgca_clean-1.cpp`, with the WELL1024a generator in `rng/` and CImg in
`CImg/` for the display build; those files are unchanged from the upstream
repository. This fork adds a browser version of the model and one study.

## LGCA Lab

[Open LGCA Lab](https://nacapule.github.io/LGCA/lgca-viz/lgca-lab.html), or open
`lgca-viz/lgca-lab.html` locally; nothing needs to be installed. It is a
JavaScript port of the C++ that reproduces it state for state: with the same seed
and settings, the parity check below finds the same lattice and the same
random-number stream. Its defaults are the C++ defaults. The alpha
normalization and the centre-plus-neighbours field used in the study are extra
options, off by default.

Parity check (needs Node.js and a C++ compiler): `node lgca-viz/verify-parity.mjs`
compiles the C++ and compares every lattice value and the full RNG state with the
JavaScript engine in 10 scenarios. See [docs/FIDELITY.md](docs/FIDELITY.md).

## Study: density × alpha (September 30)

For which alpha does the cluster density keep growing as the number of particles
grows? The alignment field is divided by M^alpha, where M is the number of
particles contributing to it: alpha 0 is the plain sum, alpha 1 the average.
Setting: polar bosons, field from the centre site plus its six neighbours,
sensitivity 6. Definitions are in [docs/PHYSICS.md](docs/PHYSICS.md).

- [Cluster density vs time at alpha 0.8 / 0.85 / 0.9 / 1](results/density-alpha-sweep-2026-09-30/figures/focus-transition-t20000-prelim.png):
  one panel per alpha, one curve per starting density (12 densities, 90×90 box).
- The same curves for [alpha 0 to 1](results/density-alpha-sweep-2026-09-30/figures/grid-time-t20000-prelim.png)
  and for [alpha 0.75 to 1 in finer steps](results/density-alpha-sweep-2026-09-30/figures/grid-time-refine-t20000-prelim.png).
- [Condensate size vs box size](results/density-alpha-sweep-2026-09-30/figures/condensate-vs-box-size-t20000.png)
  at fixed average density 2.4, box side 90 / 120 / 180 / 240. N_c is the largest
  number of particles in one channel of one site. Its late mean goes
  91 → 170 → 700 → 1887 at alpha 0.9 and 24 → 30 → 40 → 69 at alpha 1.

Curves are seed means (5 seeds; 3 for the finer alphas) up to 20,000 steps. The
numbers behind every figure are in the [study folder](results/density-alpha-sweep-2026-09-30/)
as CSV files, and `build/figures.py` redraws the figures from them.

## Resumen

LGCA Lab es una versión en JavaScript de `lgca_clean-1.cpp` que reproduce el
programa en C++ estado por estado (verificado en 10 escenarios); se abre en el
navegador sin instalar nada. El estudio del 30 de septiembre pregunta para qué alpha la densidad del cúmulo
sigue creciendo al aumentar el número de partículas. Las figuras muestran la
densidad del cúmulo contra el tiempo para 12 densidades iniciales y el tamaño del
condensado N_c contra el tamaño de la caja.
