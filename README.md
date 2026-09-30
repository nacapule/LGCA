# LGCA

Hexagonal lattice-gas cellular automaton for self-propelled particles, with
fermion (exclusion, polar alignment) and boson (unbounded occupancy) variants.
The reference code is `lgca/lgca_clean-1.cpp`, with the WELL1024a generator in
`rng/` and CImg in `CImg/` for the optional display build. Those files are
unchanged from the upstream repository.

## LGCA Lab

`lgca-viz/lgca-lab.html` is a JavaScript port of `lgca_clean-1.cpp` in a single
HTML file. With the same seed and settings it reproduces the C++ run state for
state: the same lattice at every step and the same random-number stream. To use
it, open `lgca-viz/lgca-lab.html` in a browser. Nothing needs to be installed.

The default settings are the C++ defaults. The extra options I added for the
studies below (for example the alpha normalization and the boson field that
includes the centre site) are separate toggles and do not change the default
dynamics. `docs/PHYSICS.md` describes the model and the observables,
`docs/FIDELITY.md` explains what exact agreement means and how it is checked,
and `docs/ARCHITECTURE.md` describes how the tool is built.

### Checking parity with the C++

```
cd lgca-viz && node verify-parity.mjs
```

This needs `node` and a C++ compiler (`c++`). The script compiles the actual
`lgca_clean-1.cpp`, runs 10 scenarios (fermions and bosons, several seeds,
sensitivities and densities), and compares every lattice value, the full RNG
state and the printed observables against the JavaScript engine. It ends with
"All 10 C++/JavaScript parity scenarios passed." `verify-alpha.mjs` and
`verify-neighbourhood.mjs` check the alpha and centre-site options in the same
way.

## Studies

| Date | Study | What it shows |
|---|---|---|
| Sep 10 | [Snapshot table](snapshot-table/tabla-snapshots.html) ([compare two rows](snapshot-table/compare.html)) | Coarse survey of model settings at density 0.4 on a 180×180 lattice, sensitivities 0.5 / 2 / 6, snapshots up to t = 10,000 |
| Sep 24 | [Alpha pilot](results/alpha-2026-09-24/index.html) | Normalization interpolating between SUM (alpha = 0) and AVERAGE (alpha = 1) for polar bosons; 60 runs, sizes 60 to 180, up to 10,000 steps |
| Sep 25 | [Alpha box-size refinement](results/alpha-refinement-2026-09-25/index.html) | Cluster density (mass / area) against time and box size for alpha 0.75 to 0.90; sizes 90 / 120 / 180, five seeds, 80,000 steps |
| Sep 25 | [Density vs time, with the centre site](results/density-2026-09-25/METHODS.md) | Four mean densities at alpha 0.8 / 0.9 / 1 with the centre site included in the field; the plots are in the refinement page above (Compare: mean densities) |
| Sep 30 | [Density × alpha sweep](results/density-alpha-sweep-2026-09-30/figures/) (in progress) | Twelve densities across alpha 0 to 1, with box-size runs and edge profiles; figures are updated as the runs finish |

Each results folder has the report page or figures, a `METHODS.md` where
available, and the aggregate measurements (CSV and JSON). The raw saved lattice
states from the individual runs are not included because of their size; the
settings, seeds and scripts used to produce them are in each `METHODS.md` and
`build/` folder.

## Resumen

LGCA Lab (`lgca-viz/lgca-lab.html`) es una versión en JavaScript de
`lgca_clean-1.cpp` que reproduce el programa en C++ estado por estado: con la
misma semilla da la misma red en cada paso y la misma secuencia de números
aleatorios. Se abre directamente en el navegador. La paridad se comprueba con
`cd lgca-viz && node verify-parity.mjs` (requiere node y un compilador de C++).
La tabla de arriba enlaza los estudios: la tabla de snapshots (10 sep), el
piloto de alpha (24 sep), el refinamiento por tamaño de caja (25 sep), la
densidad contra el tiempo con el sitio central (25 sep) y el barrido de
densidad × alpha (30 sep, en curso).
