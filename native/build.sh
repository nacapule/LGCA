#!/usr/bin/env bash
# native/build.sh — build the native engine's command-line program.
#
#   bash native/build.sh        ->  native/build/lgca_native
#
# It links zlib (-lz, part of macOS and of every Linux distribution) to read and write
# the sweeps' gzip-compressed lattices (`lgca_native measure` and `lgca_native run`).
#
# It also records which sources the program was built from: the SHA-256 of all the
# native/*.cpp and native/*.hpp files, joined in name order, and the flags below. The
# `run` mode writes both into each run.json's "engine" field.
#
# Compiler flags (docs/FIDELITY.md):
#   -std=c++17          the language level of the reference program
#   -O3                 full optimization; it does not change floating-point results
#   -ffp-contract=off   no fused multiply-add: clang would otherwise fuse a*b + c into
#                       one rounding, and the results would differ from the reference
#                       and from the Lab's JavaScript engine
# Never add -ffast-math or -Ofast: they allow the compiler to reorder floating-point
# arithmetic, which breaks bit-for-bit agreement with the reference.
# Set CXX to choose the compiler (default: c++). CXXFLAGS from the environment are
# ignored on purpose, so a stray global setting cannot change the arithmetic.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CXX="${CXX:-c++}"
FLAGS=(-std=c++17 -O3 -ffp-contract=off -Wall -Wextra)

# The sources' SHA-256 (sha256sum on Linux, shasum on macOS).
SOURCES=()
while IFS= read -r f; do SOURCES+=("$f"); done < <(cd "$HERE" && LC_ALL=C ls -1 -- *.cpp *.hpp | LC_ALL=C sort)
if command -v sha256sum >/dev/null 2>&1; then
    SOURCES_SHA256="$(cd "$HERE" && cat -- "${SOURCES[@]}" | sha256sum | cut -d' ' -f1)"
else
    SOURCES_SHA256="$(cd "$HERE" && cat -- "${SOURCES[@]}" | shasum -a 256 | cut -d' ' -f1)"
fi

# Built under a temporary name and then renamed over the old program, so that runs
# still using the old program are not disturbed.
mkdir -p "$HERE/build"
TMP="$HERE/build/lgca_native.tmp.$$"
trap 'rm -f "$TMP"' EXIT
"$CXX" "${FLAGS[@]}" \
    -DLGCA_SOURCES_SHA256="\"$SOURCES_SHA256\"" -DLGCA_BUILD_FLAGS="\"${FLAGS[*]}\"" \
    -o "$TMP" \
    "$HERE/lgca_cli.cpp" "$HERE/lgca_engine.cpp" \
    "$HERE/measure_cli.cpp" "$HERE/measure.cpp" "$HERE/sha256.cpp" \
    "$HERE/run_cli.cpp" "$HERE/sweep_job.cpp" "$HERE/trace_cli.cpp" \
    "$HERE/json.cpp" "$HERE/file_io.cpp" -lz
mv -f "$TMP" "$HERE/build/lgca_native"
echo "built $HERE/build/lgca_native"
