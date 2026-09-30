# How we check the JavaScript against the C++

Reference is `lgca/lgca_clean-1.cpp` with WELL1024a from `rng/`. With the same seed and
the C++ settings, the Lab keeps the same channel occupations, the same random draws in
the same order, and the same order of the arithmetic where it decides something.

From the repo root (needs Node.js and a C++ compiler):

    node lgca-viz/verify-parity.mjs

It takes the engine directly from the Lab HTML, and compiles copies of the C++ with
some prints added, in a temporary folder. In 10 cases with fermions and bosons it
compares every lattice channel and the 33 words of the random generator state at the
end; the printed observables are compared with a relative tolerance of 5e-14. It also
checks that particles are conserved, the fermion exclusion, the replay, the HTML ids and
the default values of the interface. When all is good it ends with:
"All 10 C++/JavaScript parity scenarios passed."

The 10 cases use the C++ settings. Alpha and centre + neighbours are options outside
the C++ model and are off by default. Between different computers the exact trajectory
can change if exp/log round different and this changes one decision just at the limit.
