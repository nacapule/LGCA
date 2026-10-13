// native/sha256.hpp — SHA-256 of a byte buffer, as a hex string.
//
// The sweep's run.json files identify each saved lattice by the SHA-256 of its raw bytes
// (the density and alpha sweep's JavaScript runner: createHash('sha256') over the
// Int32Array's buffer, written as 64 lowercase hex digits). The native program writes
// the same field, so it needs the same hash. This is the plain algorithm of FIPS 180-4,
// with no dependencies, checked against Node's on edge-case lengths and on saved
// lattices.

#pragma once

#include <cstddef>
#include <string>

namespace lgca {

// SHA-256 of `size` bytes at `data`, as 64 lowercase hex digits.
std::string sha256_hex(const void* data, std::size_t size);

}  // namespace lgca
