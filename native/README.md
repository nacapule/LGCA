# Lab C++ source

The Lab's simulation is implemented in C++ and compiled to WebAssembly so
it can run in a browser. The compiled engine is included in
[lgca-lab.html](../lgca-viz/lgca-lab.html); opening the Lab needs no build
tools or installation.

| File | Contents |
|---|---|
| [lgca_engine.cpp](lgca_engine.cpp) | Initialization, collisions, streaming and measurements |
| [lgca_engine.hpp](lgca_engine.hpp) | Parameters, research settings and engine interface |
| [well1024a.hpp](well1024a.hpp) | WELL1024a random generator |
| [wasm/lgca_wasm.cpp](wasm/lgca_wasm.cpp) | Functions that connect the browser to the C++ engine |
| [build-wasm.sh](build-wasm.sh) | Compiles the engine for the browser |
| [inline-wasm.mjs](inline-wasm.mjs) | Places the compiled engine in the Lab HTML |

The [model guide](../docs/PHYSICS.md) describes the collision rules and
observables. The [verification guide](../docs/FIDELITY.md) explains the
reference comparisons and numerical limits.

## Rebuild the Lab engine

For readers changing the C++ source, rebuilding needs Emscripten (`em++`
on PATH), Node.js and a C++ compiler for the reference checks.
The included engine was built with Emscripten 6.0.10; other compiler
versions can produce different module bytes.

From the repo root:

```sh
bash native/build-wasm.sh
node native/inline-wasm.mjs
node lgca-viz/verify-parity-wasm.mjs
node lgca-viz/verify-parity.mjs
```

The build produces `native/build/lgca_wasm.js`, which contains the compiled
engine and its loader. The second command embeds it in the HTML. The last
two commands check the resulting Lab against the reference model.

To compare a build with the embedded copy without changing the HTML:

```sh
node native/inline-wasm.mjs --check
```
