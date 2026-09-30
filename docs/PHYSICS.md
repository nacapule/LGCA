# PHYSICS.md — the model, exactly as implemented

The reference is `lgca/lgca_clean-1.cpp` at git HEAD (commit e216767, the advisor's
GitHub repo tailswalker/LGCA). The JS engine in `lgca-viz/lgca-lab.html` is a
bit-exact port of it. This file states the mathematics; FIDELITY.md states why the
port is exact and how that is enforced.

## Lattice and update rule

- 2-D hexagonal lattice, XDIM×YDIM sites (default 120×120), periodic boundaries.
- Each site has 6 velocity channels, one per hex direction. No rest channel.
- Channel unit vectors, index k = 0…5 (h = √3/2, computed as `sqrt(3.0)/2.0`):
  c₀=(1,0) c₁=(0.5,h) c₂=(−0.5,h) c₃=(−1,0) c₄=(−0.5,−h) c₅=(0.5,−h)
- One time step, over the whole lattice:
  1. COLLISION — each non-empty site re-picks its channel occupancies,
     conserving its particle count. Empty sites are skipped entirely.
  2. STREAMING — every particle hops one site along its channel.
- Axial coordinates: streaming destination of channel k from site (i,j) is
  DST = [(+1,0), (+1,−1), (0,−1), (−1,0), (−1,+1), (0,+1)] (k = 0…5), wrapped
  mod (XDIM, YDIM). Decreasing j corresponds to +y. The interaction
  neighbourhood (a site's 6 nearest neighbours) is
  NB = [(+1,0), (0,+1), (−1,0), (0,−1), (+1,−1), (−1,+1)].
- Kernels: J[i][j] = cᵢ·cⱼ (polar), J2[i][j] = (cᵢ·cⱼ)² (nematic).

## Two particle statistics (compile-time in C++, runtime toggle in the tool)

### FERMION — exclusion, polar alignment, Metropolis
- At most one particle per channel (site state = 6 bits).
- Collision at (x,y): field from the 6 NB neighbours only:
  m_j = Σ_{neighbours} n_j,  h_i = Σ_j J[i][j]·m_j,  E(conf) = sens·Σ_{i∈conf} h_i.
- Metropolis with 12 proposals: each proposal draws a uniformly random
  permutation of the 6 channels (Fisher–Yates), applies it to the CURRENT
  accepted configuration, and accepts iff log(u) < E_trial − E with u ~ U[0,1).
  Particle number is conserved by construction (permutations).
- dens ∈ [0,1) is the per-channel fill probability. dens ≥ 1 is outside the
  model's domain (the C++ documents this; the UI clamps typed values to 0.95
  on Reset — the engine itself reproduces the C++'s literal behaviour if driven
  directly, see FIDELITY.md).

### BOSON — unbounded occupancy, nematic, Boltzmann redraw
- Any number of particles per channel (int per channel).
- Collision at (x,y): every particle at the site independently redraws its
  channel from a softmax built from the SITE'S OWN occupancies (no neighbours):
  w_c = sens·Σ_k n_k·J2[k][c];  P(c) ∝ exp(w_c − max_c w_c).
  Sampling: r ~ U[0,1), first c with r ≤ CDF[c], fallback c = 5.
- dens ≥ 0: integer part fills every channel uniformly; fractional part is one
  extra particle per channel with probability frac(dens).

## Observables (all 0→1; N = total particles)

- polar  P = |Σ n_k c_k| / N — length of the mean velocity. Opposite
  directions cancel; 1 ⇔ common direction. Fermion order parameter.
- nematic Q = |Σ n_k (cos 2θ_k, sin 2θ_k)| / N — doubled angles, so θ and
  θ+180° reinforce: measures a common AXIS. Computed via
  cos 2θ = 2cx²−1, sin 2θ = 2·cx·cy. A perfect flock scores 1 on both P and Q.
- spatial = 1 − S/ln(XDIM·YDIM), S = −Σ_{occupied sites} p ln p, p = n_site/N.
  Even spread → 0; all particles on one site → 1. Direction-blind.
- band: for each occupied site, consider the 24 sites on rings at axial
  distances k = 1…4 in the 6 axial directions; sum doubled-angle vectors of the
  ring-site populations keyed by the RING DIRECTION index d (displacement, not
  velocity): x₃ += n·(2c_d,x²−1), y₃ += n·(2c_d,x c_d,y). Per-site score
  nsite·√(x₃²+y₃²)/(4·nband), all averaged over N. High ⇔ the mass around
  sites lies on a common axis, i.e. stripes.
  Two quirks inherited from the reference ON PURPOSE (do not "fix"):
  channel 0 is excluded from both ring counts and site weight (m = 1…5), and
  nsite accumulates the site's channels-1–5 count once per ring (×4, cancelled
  by the 4 in the denominator).
- The C++ prints, at the last step, averaged over realizations:
  sens·dens, polar±se, nematic±se, spatial(=1−ē/lnA)±se, band±se, 1/occupied±se
  with se = √(|mean² − mean-of-squares|/iters). `runBatch` reproduces this line.

## Physics toggles (defaults = lgca_clean-1.cpp exactly)

| toggle | default (C++) | alternative | provenance of the alternative |
|---|---|---|---|
| kernel | sum | average: divide the field by the particle count that built it | fermion: the commented-out `prob/=parts` in lgca_noib.cpp (~line 161). boson: ACTIVE in lgca_noib.cpp (`probi[i]/=numparts`, ~line 218, applied after max-subtraction, before exp — same order in the tool) |
| boson field from | own site | 6 NB neighbours | active in lgca_noib.cpp boson code (site term commented out) |
| boson field from | own site | own site + 6 NB neighbours (`both`) | Manik's clarification, 2026-09-25; optional seven-site research variant |
| boson alignment | nematic (J²) | polar (J) | active in lgca_noib.cpp boson code |

avg + neighbours + polar together = the ACTIVE boson collision rule of
lgca_noib.cpp (its different RNG initialization is NOT reproduced).
Note: dividing by numparts after max-subtraction equals softmax of w/numparts
(the shifted constant cancels in normalization).

## Parameter landscape (measured with the exact engine, 120², t=250, 2 seeds)

Fermion, dens 0.2 (order parameter: polar):
sens   0.2    0.3    0.35   0.4    0.45   0.6    0.8
polar  0.003  0.005  0.016  0.222  0.326  0.490  0.640
→ transition onset ≈ sens 0.4. Presets: noise 0.2 / onset 0.4 / flocking 0.8 /
dilute 1.5@d0.05.

Boson, dens 0.4 (stripes show in band + spatial, NOT in global nematic — with
the site-local kernel different stripes pick different axes):
sens     0.5    0.8    1.0    1.2    2.0
band     0.111  0.685  0.635  0.595  0.542
spatial  0.026  0.224  0.214  0.206  0.198
→ clustering onset between 0.5 and 0.8. Preset "tipping point" 1.0 is a round
marker just above onset, NOT a measured critical value. Presets: noise 0.5 /
tipping point 1.0 / stripes 2.0 / dilute 2.0@d0.1.

Saturation (zero-temperature limit): trajectories become sens-independent —
boson from sens ≈ 16, fermion from sens ≈ 32 (every accept/reject and every
softmax draw saturates; draw COUNTS are sens-independent so streams stay
aligned). The sens slider caps at 32 for this reason; the number box is
unbounded. Verified stable and conservation-exact up to sens = 10⁶.

True critical points would need longer runs, many seeds, and finite-size
scaling across lattice sizes — that is what runBatch sweeps are for.

## Power normalization experiment (2026-09-24)

Optional `kernel:"power", alpha` divides the field by M^alpha, where M is
exactly the same contributing **particle count** used by `avg`. The UI and
engine currently accept 0 ≤ alpha ≤ 2. Alpha=0 routes through the unchanged
sum arithmetic; alpha=1 through the unchanged average arithmetic. This keeps
both endpoint trajectories and WELL streams exact, including fermion and all
existing boson field/alignment combinations. Other alpha values are a new
research extension, motivated by the 23 September meeting with Manik
and Fernando. Zoom's literal `1/alpha` is inconsistent with its alpha=0 endpoint;
M^alpha is the documented working interpretation, awaiting their confirmation.
The C++ reference and all shipped physics defaults remain unchanged.

For bosons, w_c = sens * sum_k(m_k J_kc) and
p_c = exp((w_c-max(w))/M^alpha) / sum_d exp((w_d-max(w))/M^alpha).
When M=0 the field is zero and p_c=1/6. In fermions the deterministic field h
is divided by M^alpha before the existing Metropolis energy calculation.
There are no extra random draws. At fixed directional composition, the boson
score scales as sens*M^(1-alpha). This is a local identity, not a proof of
condensation or a global equivalence to a constant sensitivity rescaling.

`site` includes the entire population of the central site (including the focal
particle); `neigh` is the **six-site ring, excluding the centre**. It is not a
seven-site radius-1 disk. An isolated boson stack has zero neighbour field and
therefore redraws uniformly with `neigh` for every alpha. Do not apply the
own-site self-trapping calculation to that configuration.

As of 2026-09-25, `both` adds the centre to the six-site ring. Its channel
population is `m_k = n_centre,k + sum_neighbours n_x,k`, and the denominator
uses `M = sum_k m_k` over all seven sites, including the focal particle.
The centre has the same per-particle weight as each neighbour. This option
applies only to bosons; fermions keep their reference neighbour rule. Default
`site` and existing `neigh` arithmetic and draw order remain unchanged.
For an isolated occupied site, `both` reduces exactly to the `site` collision;
it can therefore retain an aligned moving pile that sees no field under `neigh`.
This local consequence does not establish a condensation boundary.

Manik's clarified diagnostic (2026-09-25): compare cluster density against time
for successively larger total populations, seek relaxed plateaus, then compare
their heights. Current study increases N by enlarging the box at fixed mean
density; a fixed-box density sweep is a separate axis. Rising plateau heights
over finite sizes are evidence to investigate, not a proof of divergence or a
measured critical alpha. All completed September24/25 sweeps predate `both`.

Added display diagnostics: n_max=max_x sum_k n_xk; peak fraction=n_max/N;
A_eff=N²/sum_x n_x² (in sites). A_eff is a participation area of the entire
population, not the geometrical area of the largest cluster. For zero particles
these diagnostic fractions are displayed as 0; original C++ observables retain
their original convention. The log density view maps log(1+n_x)/log(1+ceiling)
to a fixed palette and clips at the ceiling; it never rescales by the instantaneous
maximum. `dens` remains an initialization parameter per channel; actual mean
site occupancy is N/(WH), approximately 6*dens.

`verify-alpha.mjs` checks legacy/alpha endpoint state and RNG equality,
intermediate-alpha invariants, an independently evaluated asymmetric-field CDF,
zero sensitivity and empty-neighbour fields; plus 24 checkpointed comparisons
with real C++ on a 96×60 lattice over four 300-step stress runs.
`verify-neighbourhood.mjs` additionally checks the seven-site field and its
normalization against independent CDF/sample calculations at periodic corners,
the isolated-site limit, conservation, RNG draw counts, endpoint equivalence,
replay, and unchanged old variants against the pinned study engine.
