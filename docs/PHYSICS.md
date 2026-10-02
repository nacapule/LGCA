# Model and observables

The model lives on a hexagonal lattice with six directions at each site.
At each whole time step, particles first choose their directions and then
move one site. There is no rest direction. Boundaries are periodic: a
particle leaving one edge enters through the opposite edge.

The input `dens` is **particles per channel**, not per site. There are six
channels, so the expected starting density is `6 × dens` particles per
site. The actual mean density is `rho = N / (W × H)`, where N is the
conserved particle total.

## Fermions and bosons

| Model | Occupancy | Reference collision rule |
|---|---|---|
| Fermion | At most one particle per channel | Polar alignment with the six neighbouring sites; 12 Metropolis proposals that permute the occupied channels |
| Boson | Multiple particles per channel | Each particle independently draws a direction from the site's own nematic alignment field |

For fermions, `dens` is a fill probability in `[0, 1)`.
For bosons, each channel starts with `floor(dens)` particles, plus one more
with probability `dens − floor(dens)`, for nonnegative density.

Channel k has unit vector `c_k = (cos(kπ/3), sin(kπ/3))`, for k = 0…5.
Polar alignment uses `J_kc = c_k · c_c`: opposite directions cancel.
Nematic alignment uses `J_kc = (c_k · c_c)²`: opposite directions reinforce
the same axis.

A fermion proposal is accepted if `log(u) < E_trial − E`,
where `E = sensitivity × sum of the alignment field over occupied channels`
and `u` is a WELL1024a draw in `[0, 1)`.

## Field settings and alpha

A field counts the particles whose directions influence a collision.

| Boson field | Sites counted |
|---|---|
| Own site | The central site, including the particle choosing a direction |
| Neighbours only | The six surrounding sites, excluding the centre |
| Site plus neighbours | All seven sites, with equal per-particle weights |

Fermions always use the six neighbours. Boson field and alignment choices
are optional; the reference defaults are **own site, nematic, sum**.

Let `m_k` be the particles in channel k across the sites contributing to
the field, and `M = Σ_k m_k`. A boson's directional score is

```text
w_c = sensitivity × Σ_k m_k J_kc
```

With power normalization, its probability to choose channel c is

```text
p_c = exp((w_c − max_d w_d) / M^alpha)
      / Σ_d exp((w_d − max_e w_e) / M^alpha)
```

| Kernel | Division | Meaning |
|---|---|---|
| Sum | None | More contributing particles make a stronger field |
| Average | M | Field per contributing particle |
| Power | M^alpha | Interpolation: alpha 0 gives sum; alpha 1 gives average |

At fixed directional proportions, the normalized score scales as
`sensitivity × M^(1−alpha)`. This describes the local field; it does not
by itself prove condensation.

If M is zero, all six boson directions have equal probability.
For fermions, the neighbour field is divided by the same normalization
before the Metropolis test.

The study uses polar bosons, site plus neighbours, sensitivity 6, and power
normalization. Its rules differ from the reference defaults.

## Measurements in the Lab

| Observable | Definition |
|---|---|
| Polar | `|Σ n_k c_k| / N`. Close to 1 when particles share a direction. |
| Nematic | `|Σ n_k (cos 2θ_k, sin 2θ_k)| / N`. Close to 1 when particles share an axis, including opposite directions. |
| Spatial | `1 − S / ln(W×H)`, with `S = −Σ p_x ln p_x` and `p_x = n_x/N`. Zero for an even spread; one for all particles at one site. |
| Band | Measures whether surrounding mass lies along an axis. It samples distances 1–4 along six axial directions around each occupied site. |

Band counts **all six channels**. The reference omits channel 0.
Older band values use that earlier definition; the change affects only
this measurement.

## Measurements in the study

Let `n_x = Σ_k n_xk` be the population of site x.

A site belongs to the cluster search if
`n_x ≥ max(1, ceil(2 × rho))`. Adjacent qualifying sites form a component,
using all six neighbours and periodic boundaries. We select the component
with the greatest particle mass.

| Measurement | Definition |
|---|---|
| Cluster mass | Particles in the selected component |
| Cluster area | Number of sites in that component |
| Cluster density | Mass / area, in particles per site; zero if no component qualifies |
| Condensate occupation N_c | `max over x,k of n_xk`: the largest occupation of a single channel at a single site |

N_c is different from the total population of the most occupied site and
from the mass or density of a connected cluster.

The published late values average steps 15,100–20,000 within each run,
then average those run means across seeds. Ranges show the smallest and
largest seed means. A flat part of a finite curve is not enough to establish
a limiting density.
