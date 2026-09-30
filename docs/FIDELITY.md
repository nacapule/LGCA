# FIDELITY.md — the bit-exactness contract and how not to break it

The central requirement of this project: the tool's default-settings dynamics must
be indistinguishable from `lgca/lgca_clean-1.cpp` — not statistically, but
STATE-FOR-STATE. Same seed + same settings ⇒ the identical sequence of lattice
states and the identical WELL1024a stream as the corresponding C++ build
(MODEL=FERMION or BOSON). Any divergence from the reference mathematics is the
project's catastrophic failure mode. When fidelity and anything else conflict,
fidelity wins.

## What the contract rests on (all verified 2026-08-11)

1. RNG: WELL1024a ported constant-for-constant from `rng/WELL1024a.c`
   (R=32 words, M1=3, M2=24, M3=10; shifts 8/19/14 and 11/7/13;
   FACT = 2.32830643653869628906e-10; index walk state_i=(state_i+31)&31).
   Seeding: init[i]=seed; seed = seed·1664525 + 1013904223 (32-bit unsigned,
   via Math.imul + >>>0). `makeWell` also exposes getState/setState (33 words:
   32 state + index) for replay — a pure addition, no effect on the stream.
2. Draw ORDER and draw COUNT, per step, identical to the C++:
   - init: one draw per (i,j,k), i outer, j inner, k innermost — even when
     frac = 0.
   - collision loop: i outer, j inner; EMPTY SITES CONSUME NO DRAWS.
   - fermion, per non-empty site: 12 proposals × (5 Fisher–Yates draws
     `(rng()*(k+1))|0` for k=5…1 + 1 acceptance draw) = 72 draws, always,
     regardless of outcomes.
   - boson, per non-empty site: CDF is deterministic; exactly n draws
     (n = site population).
   - streaming, measurement, rendering, UI: ZERO draws. The seed-dice button
     uses Math.random — never sim.rng.
3. Floating-point ARITHMETIC ORDER in everything that feeds a decision or a
   published number:
   - J tables built from `Math.sqrt(3.0)/2.0` (sqrt is IEEE-correctly-rounded,
     so bits match the C++'s std::sqrt exactly). Literal tables like 0.25/0.5
     would be 1 ulp off in places and DO diverge — this was hit in practice.
   - field/energy sums in the same nesting order as the C++ loops.
   - measure(): i-outer/j-inner site order; √(x²+y²) via Math.sqrt of the sum
     (the C++'s sqrt(sqr+sqr)) — NOT Math.hypot (differently rounded).
   - stderr: mean*mean, not mean**2 (** may lower to a pow call).
   - boson softmax: subtract max, then (avg kernel only) divide, then exp —
     the exact sequence of the corresponding C++.
4. Integer state: channel occupancies are Int32Array, matching C++ `int`.
   (Uint16 would wrap silently at 65,536 in a clustered boson channel.)
5. The C++ side quirks are reproduced, not repaired: band's channel-0
   exclusion and ×4 site weight (PHYSICS.md); fermion dens ≥ 1 garbage-in
   (nuparts inflated by base·6·X·Y while the lattice starts empty — the UI
   Reset clamps typed input to 0.95, the engine drivable past it on purpose);
   dens = 0 → NaN observables (0/0), as in the C++.

## Rules for any future change

- NEVER add, remove, or reorder RNG draws in the simulation path of the
  C++-default modes. New stochastic features must draw from a SEPARATE
  generator (as tracers once did) or exist only behind non-default toggles.
- NEVER reorder floating-point accumulation in decision paths or observables
  of the default modes. Optimizations must produce the identical operation
  sequence (bit-packing, hoisting invariant INTEGER work: fine; re-associating
  float sums: not fine).
- New physics variants: labelled toggle, default = C++ behaviour, provenance
  documented in PHYSICS.md. The parity suite asserts the shipped defaults.
- The engine lives in the `<script id="engineSrc">` block; the parity script
  extracts THAT block, so the shipped engine is what gets tested. Keep the id.
- Keep occ/src/tmp Int32Array (the suite asserts it).

## Verification protocol

After ANY edit to `lgca-viz/lgca-lab.html` — engine OR ui:

    cd lgca-viz && node verify-parity.mjs

must end with "All 10 C++/JavaScript parity scenarios passed." It compiles the
real C++ (needs `node` ≥ 18 and `c++`), then checks: full-lattice integer
equality AND all-33-word RNG-state equality per scenario, metrics to 5e-14,
per-step conservation + fermion exclusion, replay round-trip exactness,
HTML/script sanity (parses, unique ids, referenced ids exist, C++-default
toggle/values unchanged). It catches UI regressions too (a deleted element id,
a changed default) — run it even for "cosmetic" changes.

For ENGINE-level changes, additionally run longer adversarial checks in the
style of the 2026-08-11 audit: checkpointed state+RNG dumps every ~50 steps
against an instrumented C++ build, several hundred steps, high-sens/high-dens
stress, and a RECTANGULAR lattice (patch XDIM≠YDIM, e.g. 96×60 — square
lattices cannot expose axis swaps). Method: prepend WELL1024a.c to the .cpp
(strip its `#include ../rng/`, #undef its macros), patch XDIM/YDIM/MODEL,
insert a dump block in the time loop, compile `c++ -std=c++17 -O2`.
Note the makefile's include paths point at $(HOME)/Documents/lgca/ and do not
work in this clone; the parity script's self-contained build is the way.

## The one bounded caveat: exp/log rounding

Everything except Math.exp and Math.log is bit-identical across platforms.
Those two differ between math libraries in the last bit. Measured V8 (node 26,
arm64) vs Apple libm on the live argument domains: log differs on 6.6% of
draws, exp on 9.8%, ALWAYS by exactly 1 ulp. A trajectory can only diverge if
a decision threshold lands inside that 1-ulp gap: probability ~1e-16 per
decision; zero flips observed in ~1e8 live draws. Consequences:
- On this machine, V8-vs-C++ parity holds in practice for full runs.
- Other JS engines (Safari/JSC, Firefox) have their own exp/log: identical
  stochastic dynamics, same 1-ulp analysis, but bitwise trajectories may
  diverge from the C++ at some point in a long run. The physics is unaffected.
- Gold-standard option if universal bitwise reproducibility is ever wanted:
  adopt a correctly-rounded exp/log (e.g. CORE-MATH) in BOTH the C++ and the
  JS. Then every platform agrees forever.

## Audit record

2026-08-11 — full independent audit:
- Line-by-line map of the whole engine against the pristine HEAD C++; RNG
  constants against WELL1024a.c; toggle provenance against lgca_noib.cpp.
- Shipped suite 10/10. Independent harness: 6 scenarios, 1700 steps total,
  34 checkpoints (every lattice integer + all RNG words), including 500-step
  runs in both physical regimes, sens 8 / dens 2.25 and dens 0.95 stress, and
  96×60 rectangular — ALL bit-exact; metrics ≤ 1.1e-16 relative.
- libm quantification as above. In-browser checks: live-UI time travel
  round-trip bit-exact; conservation and exclusion hold at every step.
- Sensitivity range extension verified safe to sens = 10⁶ (conservation,
  exclusion, finiteness; saturation from ≈16/≈32 — see PHYSICS.md).
