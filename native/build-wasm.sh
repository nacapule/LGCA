#!/usr/bin/env bash
# native/build-wasm.sh — build the engine as WebAssembly for the Lab and for Node.
#
#   bash native/build-wasm.sh   ->  native/build/lgca_wasm.js
#
# The output is ONE JavaScript file (the "glue"): the compiled WebAssembly is embedded
# in it as a base64 string, so it can be inlined into lgca-viz/lgca-lab.html as a text
# block and still work from a file:// page, inside a Blob-URL worker, or in Node.
#
# It compiles the same engine the native program runs (native/lgca_engine.cpp, unchanged)
# plus a plain C interface for JavaScript (native/wasm/lgca_wasm.cpp, documented there).
# Needs Emscripten's em++ (`brew install emscripten`; em++ is emcc for C++, which also
# links the C++ runtime).
#
# Compiler flags (docs/FIDELITY.md), the same as native/build.sh:
#   -std=c++17          the language level of the reference program
#   -O3                 full optimization; it does not change floating-point results
#   -ffp-contract=off   no fused multiply-add (WebAssembly has none in its base
#                       instruction set, but the flag keeps the rule explicit)
# Never add -ffast-math or -Ofast: they allow the compiler to reorder floating-point
# arithmetic, which breaks bit-for-bit agreement with the reference.
#
# Emscripten settings:
#   MODULARIZE=1, EXPORT_NAME=createLgcaWasm
#       the glue defines one factory function, createLgcaWasm(), returning a Promise of
#       a fresh module instance. A classic script, NOT -sEXPORT_ES6: the ES-module glue
#       uses import.meta, which does not load inside a Blob-URL worker or from an
#       inlined text block. Node loads the classic glue with require().
#   ENVIRONMENT=web,worker,node
#       the places the glue must run: the Lab page, its worker, and the Node checks.
#   ALLOW_MEMORY_GROWTH=1
#       the module's memory grows to fit the lattice (up to 2 GiB). A JavaScript view
#       on that memory (HEAP32 and the others) must be re-read after any call that
#       may allocate, because growing replaces the underlying buffer.
#   SINGLE_FILE=1, SINGLE_FILE_BINARY_ENCODE=0
#       embed the WebAssembly in the glue as base64. Emscripten 6 otherwise uses a
#       raw UTF-8 embedding that breaks if a page is read in any other character set.
#   FILESYSTEM=0
#       no emulated file system: the engine reads and writes no files.
#   EXPORTED_FUNCTIONS, EXPORTED_RUNTIME_METHODS
#       the C functions of lgca_wasm.cpp, malloc/free for passing arrays, and the
#       memory views and string reader the JavaScript side uses.
# Set EMXX to choose the compiler (default: em++). CXXFLAGS and EMCC_CFLAGS from the
# environment are ignored on purpose, so a stray global setting cannot change the
# arithmetic.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EMXX="${EMXX:-em++}"
unset EMCC_CFLAGS

if ! command -v "$EMXX" >/dev/null 2>&1; then
    echo "em++ not found: install Emscripten (brew install emscripten)" >&2
    exit 1
fi

FLAGS=(-std=c++17 -O3 -ffp-contract=off -Wall -Wextra)

EXPORTS=(
    _malloc _free
    _lgca_create _lgca_destroy _lgca_error
    _lgca_init_lattice _lgca_step
    _lgca_occ_ptr _lgca_src_ptr _lgca_cell_count _lgca_set_occ
    _lgca_get_rng _lgca_set_rng
    _lgca_measure _lgca_band _lgca_nuparts _lgca_t _lgca_set_nuparts _lgca_set_t
    _lgca_set_options
)
RUNTIME=(HEAP32 HEAPU32 HEAPF64 HEAPU8 UTF8ToString)

join() { local IFS=,; echo "$*"; }

mkdir -p "$HERE/build"
"$EMXX" "${FLAGS[@]}" \
    -sMODULARIZE=1 -sEXPORT_NAME=createLgcaWasm \
    -sENVIRONMENT=web,worker,node \
    -sALLOW_MEMORY_GROWTH=1 \
    -sSINGLE_FILE=1 -sSINGLE_FILE_BINARY_ENCODE=0 \
    -sFILESYSTEM=0 \
    -sEXPORTED_FUNCTIONS="$(join "${EXPORTS[@]}")" \
    -sEXPORTED_RUNTIME_METHODS="$(join "${RUNTIME[@]}")" \
    -o "$HERE/build/lgca_wasm.js" \
    "$HERE/wasm/lgca_wasm.cpp" "$HERE/lgca_engine.cpp"

BYTES="$(wc -c < "$HERE/build/lgca_wasm.js" | tr -d ' ')"
echo "built $HERE/build/lgca_wasm.js ($BYTES bytes, $("$EMXX" --version | head -n 1))"
