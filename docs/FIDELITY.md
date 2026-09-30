# Checking the JavaScript port

The reference is `lgca/lgca_clean-1.cpp`, with WELL1024a from `rng/`.
For matching seeds and C++ settings, the Lab preserves channel occupancies,
random draws and their order, and the arithmetic order used in decisions.

Run from the repository root (requires Node.js and a C++ compiler):

    node lgca-viz/verify-parity.mjs

The checker extracts the engine directly from the shipped Lab HTML and compiles
instrumented copies of the reference C++ in a temporary directory.
Ten fermion and boson scenarios compare every lattice channel and all 33 RNG
state words at their endpoints; printed observables use a 5e-14 scaled tolerance.
It also checks conservation, fermion exclusion, replay, HTML IDs and UI defaults.
Success ends with: "All 10 C++/JavaScript parity scenarios passed."

The parity scenarios use the C++ settings. The alpha and centre-plus-neighbours
options go beyond the C++ model and are off by default.
Exact cross-platform trajectories can differ when exp/log rounding changes
a decision at its threshold.
