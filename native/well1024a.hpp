// native/well1024a.hpp — the WELL1024a random number generator as a small struct.
//
// A constant-for-constant port of rng/WELL1024a.c, the generator the reference program
// (lgca/lgca_clean-1.cpp) links. The original keeps its state in file-level globals; here
// the same 32 state words and the same index live in a struct, so one program can hold
// several independent generators (tests, sweeps, later WebAssembly).
//
// WELL1024a: F. Panneton, P. L'Ecuyer and M. Matsumoto, "Improved long-period generators
// based on linear recurrences modulo 2", ACM Transactions on Mathematical Software 32
// (2006). The original code (rng/WELL1024a.c) is by Francois Panneton and Pierre L'Ecuyer
// (University of Montreal) and Makoto Matsumoto (Hiroshima University); its notice
// allows free personal, academic and non-commercial use.
//
// What must not change (docs/FIDELITY.md): the constants R = 32, M1 = 3, M2 = 24,
// M3 = 10, the shifts 8 / 19 / 14 and 11 / 7 / 13, the factor FACT, the index walk
// index = (index + 31) & 31, and the seeding rule of the reference main(). The first
// 10^6 doubles are checked against the original C code.

#pragma once

#include <cstdint>

namespace lgca {

struct Well1024a {
    static constexpr int R  = 32;   // number of 32-bit state words
    static constexpr int M1 = 3;
    static constexpr int M2 = 24;
    static constexpr int M3 = 10;
    // 2^-32: turns a 32-bit word into a double in [0, 1).
    static constexpr double FACT = 2.32830643653869628906e-10;

    uint32_t state[R] = {};   // STATE[] in the original
    uint32_t index = 0;       // state_i in the original

    // InitWELLRNG1024a(init): copy the 32 words, start at index 0.
    void init(const uint32_t init_words[R]) {
        index = 0;
        for (int j = 0; j < R; ++j) state[j] = init_words[j];
    }

    // The seeding of the reference main(): all 32 words from one 32-bit seed through a
    // one-line linear congruential generator (unsigned arithmetic wraps modulo 2^32).
    void seed_like_reference(uint32_t seed) {
        uint32_t init_words[R];
        for (int i = 0; i < R; ++i) {
            init_words[i] = seed;
            seed = seed * 1664525u + 1013904223u;
        }
        init(init_words);
    }

    // WELLRNG1024a(): advance the generator and return a double in [0, 1).
    // The original's macros, written out:
    //   MAT0POS(t, v) = v ^ (v >> t)        MAT0NEG(t, v) = v ^ (v << -t)
    //   V0 = state[i], VM1 = state[(i+M1)&31], VM2 = state[(i+M2)&31],
    //   VM3 = state[(i+M3)&31], VRm1 = newV0 = state[(i+31)&31], newV1 = state[i].
    double next() {
        const uint32_t v0   = state[index];
        const uint32_t vm1  = state[(index + M1) & 0x1fu];
        const uint32_t vm2  = state[(index + M2) & 0x1fu];
        const uint32_t vm3  = state[(index + M3) & 0x1fu];
        const uint32_t vrm1 = state[(index + 31) & 0x1fu];

        const uint32_t z0 = vrm1;
        const uint32_t z1 = v0 ^ (vm1 ^ (vm1 >> 8));                   // Identity(V0) ^ MAT0POS(8, VM1)
        const uint32_t z2 = (vm2 ^ (vm2 << 19)) ^ (vm3 ^ (vm3 << 14)); // MAT0NEG(-19, VM2) ^ MAT0NEG(-14, VM3)

        state[index] = z1 ^ z2;                                        // newV1
        state[(index + 31) & 0x1fu] = (z0 ^ (z0 << 11))                // newV0 = MAT0NEG(-11, z0)
                                    ^ (z1 ^ (z1 << 7))                 //       ^ MAT0NEG(-7, z1)
                                    ^ (z2 ^ (z2 << 13));               //       ^ MAT0NEG(-13, z2)
        index = (index + 31) & 0x1fu;
        return static_cast<double>(state[index]) * FACT;
    }

    // The full generator state in the layout of the Lab's JavaScript getState():
    // the 32 state words, then the index (33 words).
    void get_state(uint32_t out[R + 1]) const {
        for (int j = 0; j < R; ++j) out[j] = state[j];
        out[R] = index;
    }

    // Restore a state saved by get_state() (or by the Lab's getState()). A saved index is
    // always below 32; the mask only keeps a corrupt input from indexing out of bounds.
    void set_state(const uint32_t in[R + 1]) {
        for (int j = 0; j < R; ++j) state[j] = in[j];
        index = in[R] & 0x1fu;
    }
};

}  // namespace lgca
