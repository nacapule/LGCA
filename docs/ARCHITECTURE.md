# ARCHITECTURE.md — repo layout and how the tool is built

## Repository

Fork of the advisor's https://github.com/tailswalker/LGCA (upstream commit
e216767). The upstream C++, RNG and CImg files are unchanged; everything else
listed below was added on top.

```
LGCA/
├── docs/                      PHYSICS.md · FIDELITY.md · ARCHITECTURE.md
├── lgca/
│   ├── lgca_clean-1.cpp       THE REFERENCE (see PHYSICS.md)
│   ├── lgca_noib.cpp          older variant; source of the physics toggles
│   └── makefile               advisor's; include paths don't fit this clone
├── rng/WELL1024a.c/.h         reference RNG
├── CImg/                      graphics lib for the C++'s optional -DCIMG build
├── lgca-viz/
│   ├── lgca-lab.html          THE TOOL — one self-contained file, no deps
│   └── verify-parity.mjs      C++↔JS parity suite (see FIDELITY.md)
├── results/                   dated studies (report pages, figures, aggregates)
└── snapshot-table/            snapshot table of the parameter landscape (Sep 10)
```

## The tool (`lgca-lab.html`) — one file, three layers

Hard constraint: stays a single self-contained HTML file. No build step, no
frameworks, no external requests. Open it in a browser and it runs.

### Engine — `<script id="engineSrc">` (the parity-tested block)
- `makeWell(seed)` — WELL1024a; gen() → double in [0,1); getState/setState
  (33-word Uint32Array) for replay.
- `class LGCA` — constructor takes {model, W, H, dens, seed, kernel,
  bosonField, bosonAlign}. State: `occ` (post-streaming, the current lattice),
  `src` (post-collision, pre-streaming — what the glide animation renders),
  `tmp` (scratch), all Int32Array of W·H·6, indexed `(j*W+i)*6+k`.
  Methods: initLattice, neighbourField, collideFermion, collideBoson,
  step(sens), measure(), band(). Loop orders and draw discipline: FIDELITY.md.
- `runBatch({model,tsteps,iters,sens,dens,seed,W,H,...})` — mirrors the C++
  main() realization loop on ONE rng stream; prints/returns the same output
  line. `dumpState({...})` — raw lattice dump for validation. Both on window —
  usable from the browser console for sweeps.

### UI — second `<script>` block
- P = current params; `sim` = live LGCA instance. makeSim/reset build it from
  the controls (reset clamps typed density: <0→0; fermion ≥1→0.95).
- Main loop: requestAnimationFrame; speed slider (SPEEDS array, Infinity =
  12 ms step budget per frame); `phase` ∈ [0,1] is the glide fraction.
- Rendering: particles mode draws `src` gliding by `phase` toward DST
  (interpolated motion IS the streaming semantics — see the in-tool notes);
  raw mode draws `occ` stroboscopically; fields mode draws `occ` into a 2·W×2·H
  offscreen image (density / classic C++ blue / direction hue / nematic hue)
  with optional axis/flux overlays. Sheared-torus projection: site (i,j) at
  x=(i+j/2)·s0 (mod world), y=j·(√3/2)·s0; wheel-zoom, drag-pan, 0/dblclick
  reset. Site inspector tooltip on hover (channel rose).
- Replay/time travel: ring buffer `rb` of {src copy, rng snapshot, tick},
  capacity ≈ 24 MB (min 20, max 300 slots). Slider scrubs VIEW-ONLY
  (replay.on); transport buttons MOVE the sim: jumpTo(k) restores lattice
  (re-streaming the buffered src), rng state, truncates hist and buffer —
  stepping forward again is bit-identical (verified). Param change after
  rewind = deliberate counterfactual branch.
- Chart: `CH` pad constants shared with the hover handler (do not duplicate
  numbers); fine grid — y minors 0.05/majors 0.25 labelled; x nice-number
  (1-2-5) majors with k/M compact labels, minors when ≥4 px and ≥1 step;
  canvas height 186. Series: polar/nematic/spatial each step; band every 10
  steps (checkbox), stored in hist.bt/bv.
- CSV button: exports step,polar,nematic,spatial,band (band sparse); clipboard
  first, honest fallback to a .csv download ("saved .csv") on denial.
- Cards: the right column's five sections are `<details class="card"
  data-card=…>` — collapsible (summary click) and drag-reorderable (⠿ handle);
  layout persists in localStorage key `lgcaCards` {order, open}.
- Transport row: Play | −1 | +1 | −N | [N] | +N | Reset (Reset last,
  amber-tinted via --reset-* tokens, both themes ≥5:1 text contrast).
  Keys: space play/pause · , / . step back/ahead (ahead animates) · [ / ]
  jump −N/+N · 0 reset view.
- Theming: light/dark via prefers-color-scheme variables on :root.

### Removed features
- PNG snapshot + WebM recording (2026-08-11), including --danger tokens.
- "Boson jiggle" fidelity note; "back is limited to the replay buffer" note
  (the −1/−N tooltips still carry the buffer-limit info).
- v1 (`lgca-explorer.html`, boson-only) was deleted; v2 renamed to
  lgca-lab.html. verify-v2-parity.mjs → verify-parity.mjs.

## Performance facts (measured)

- Single JS thread: sim + rendering share one core; Chrome Helper at 100% =
  one core saturated. 120²: boson ~2–3 ms/step, fermion ~16 ms (72 draws + 12
  energy evaluations per occupied site); 360² fermion ~132 ms.
- Possible future directions, in order of value: (1) Web Worker for the sim
  (responsiveness, zero fidelity cost); (2) parallel SWEEPS — one worker per
  core, each an independent exact run (the research win); (3) intra-run
  parallelism conflicts with the single sequential RNG stream (only exotic
  WELL jump-ahead could keep exactness); (4) WASM ~1.5–3× (exp/log caveat);
  GPU = loosest parity. Per-core micro-opts must preserve op order
  (FIDELITY.md); bit-packing the fermion state is the safe candidate.

## Alpha and density additions (2026-09-24)

The general single-file Lab now exposes optional `kernel:"power"` with alpha,
plus the existing sum/average modes. Alpha is serialized with recipes and
copied to the live engine; its endpoints preserve the old engine arithmetic.
`?study=alpha&alpha=...&size=...&dens=...&sens=...&seed=...&field=site|neigh|both` selects a
validated explicit polar-boson study setup and pauses. Ordinary launches retain
all reference defaults. This link starts fresh (the existing Reset advances to
step 1); it is not a late-time state loader.
`dens` is optional, defaults to .4 per channel, and must be finite and nonnegative.

Density is now directly selectable. Log density has a user-set fixed count
ceiling; the old density/blue/direction/axis fields remain available. Display-only
peak/fraction/participation-area metrics and per-step active physics fields are
included in CSV. Density snapshots retain PNG and occupancy/RNG/metadata in the
current page session and offer file exports. Saved settings remain recipe-only.
Alpha has .01 slider/numeric steps and explicit ±.1 buttons. The boson field
control includes own-site, neighbours-only, and site-plus-neighbours (`both`).
Recipes, density exports and CSV preserve that field setting.

`results/alpha-2026-09-24` is an independent analysis artifact, not an alternate
engine. The completed September24 and25 studies now load the original engine
from `results/alpha-2026-09-24/build/engine-source.js`, whose SHA-256 matches
their recorded hashes, so adding the optional `both` mode does not break their
reproduction/resume scripts. New studies should extract and record the current
Lab engine. `build/lgca-lab.before.html` is a frozen copy of the Lab from before the alpha
changes, kept as a regression reference. The report and its images/figures can
be opened directly. Cross-links to the Lab use repo paths.
