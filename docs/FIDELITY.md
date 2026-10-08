# Model verification

These checks are for changes to the simulation code. Opening and using the
Lab requires only its HTML file and a browser.

The Lab runs the C++ engine in [native/](../native/README.md), compiled to
WebAssembly and embedded in its HTML file. The reference model is
[lgca_clean-1.cpp](../lgca/lgca_clean-1.cpp), with the WELL1024a generator
in [rng/](../rng/).

For the same seed and reference settings, the check compares every lattice
channel and every word of the random generator state. It tests the compiled
engine included in the Lab; no rebuild is needed to run the check.

## Check against the reference

From the repo root, with Node.js and a C++ compiler:

```sh
node lgca-viz/verify-parity-wasm.mjs
```

The script reads the compiled engine from the Lab HTML and builds temporary
copies of the reference C++ with state-printing probes. It does not replace
the source files.

| Check | Comparison |
|---|---|
| Ten fermion and boson cases | Every lattice channel and all 33 generator words against the reference |
| Printed measurements | Error no greater than `5e-14 × max(1, \|reference value\|)` |
| Replay | Restoring the lattice and generator reproduces the same continuation |
| Research settings | Step-by-step agreement with an independent implementation across 60 combinations |
| Buffered measurements and rewinds | Worker replies agree with fresh runs from the seed |

Success ends with:

```text
All 10 C++/WebAssembly parity scenarios passed.
```

## The band correction

Our band measurement counts all six channels. The reference starts at
channel 1 and skips channel 0. The verifier applies the same one-line
correction to its temporary C++ copies.

This changes band values, not particle motion or random draws. Values made
with the old definition should be kept separate from corrected values.
The omitted channel can bias band upward or downward, depending on the state.

## Scope of verification

For the checked seeds and settings, the complete lattice and random generator
states match. This is stronger than comparing pictures or averages.

It does not prove equality for every possible parameter, seed or machine.
Different math libraries can round `exp`, `log` and `pow` differently in the
last bit. If that changes a decision near its threshold, the later trajectory
can differ. The compiled browser engine carries its own math routines, so
those routines are the same across browsers running the same module.

Alpha normalization and alternative fields are optional changes to the
reference model. Their comparison checks agreement between implementations;
it does not establish a physical result for every research setting.

## Additional checks and engine changes

The HTML also contains a JavaScript fallback for cases where the compiled
module cannot load. Its comparison check is:

```sh
node lgca-viz/verify-parity.mjs
```

It checks ten reference cases, conservation, fermion exclusion, replay and
default settings. Run both checks after any change to the Lab HTML.
Comparisons between implementations allow a rounding difference of
`1e-12 × max(1, |value|)` for entropy and spatial order; lattice and
generator states must match exactly.

For changes to the simulation:

- Preserve the WELL1024a state and the count and order of all draws.
- Preserve the order of floating-point sums and divisions.
- Keep channel populations in signed 32-bit storage.
- Compile comparisons with `-ffp-contract=off`; do not use
  `-ffast-math` or `-Ofast`.
- Change the C++ source and [rebuild the embedded engine](../native/README.md#rebuild-the-lab-engine), rather than editing the generated HTML block.
