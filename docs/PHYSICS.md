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

Here `n_xk` is the number of particles at site x in channel k,
`n_x = Σ_k n_xk`, `n_k = Σ_x n_xk`, and `θ_k = kπ/3` is channel k's angle.

| Observable | Definition |
|---|---|
| Polar | `\|Σ n_k c_k\| / N`. Close to 1 when particles share a direction. |
| Nematic | `\|Σ n_k (cos 2θ_k, sin 2θ_k)\| / N`. Close to 1 when particles share an axis, including opposite directions. |
| Spatial | `1 − S / ln(W×H)`, with `S = −Σ p_x ln p_x` and `p_x = n_x/N`. Zero for an even spread; one for all particles at one site. |
| Band | Whether the mass around each occupied site lies along an axis, weighted by the site's population (below) |

For each occupied site x, the band samples the sites `y = x + r·e_d` at
distances r = 1…4 along each of the six lattice directions d (24 sites),
where `e_d` is the hop of channel d:

```text
v_x = Σ_{r,d} n_y (cos 2θ_d, sin 2θ_d)
m_x = Σ_{r,d} n_y
b_x = |v_x| / m_x        (0 when m_x = 0)
Band = Σ_x n_x b_x / N
```

Doubling the angle makes opposite directions count as the same axis. On a
small periodic lattice two offsets can reach the same site; each is counted.

Band counts **all six channels**. The reference omits channel 0.
Older band values use that earlier definition; the change affects only
this measurement.

## Measurements in the study

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

In the density and alpha study, late values average steps 15,100–20,000
within each run, then average those run means across seeds. Ranges show
the smallest and largest seed means. A flat part of a finite curve is not
enough to establish a limiting density.

## Islands

The island study uses connected
occupied sites with no density cutoff. Two sites are connected when both
hold at least one particle and they are neighbours (six neighbours,
periodic boundaries). An **island** is a group of connected occupied sites
together with every occupied site connected to it. Its density is its
particle count divided by its number of sites.

| Island | Definition |
|---|---|
| Busiest site's island | The island that contains the site with the most particles. A tie keeps the first site in index order. |
| Heaviest island | The island with the most particles |
| Densest island | The island with the most particles per site. It can be a single site. |

The cluster above keeps only sites at or above the cutoff, while an island
keeps every occupied site it connects. So the cluster can be denser than
any island.

At 2.4 particles per site a new lattice has about 95% of its sites
occupied, more than the site-percolation threshold of the six-neighbour
lattice, which is 1/2. In the study's runs all of these sites form one
island at step 0. The busiest site's island separates once the gas between
piles thins out.

Late island values average each run's samples at steps `0.75 T < t ≤ T`,
then average those run means over the seeds. T is the last step of the
runs; the larger-box runs also have late values at T = 20,000. The densest
island is measured only on saved lattices, so its late value uses the
lattice saved at step T.

### Island CSV columns

The island study's three CSV files use these columns.

| Column | Meaning |
|---|---|
| `box` | `denser`: 90 × 90 box, N set by `dens`. `larger`: `dens` 0.4, N set by the box side |
| `alpha`, `side`, `dens` | Alpha, box side L (L × L sites), starting `dens` per channel |
| `N` | Particle count of a saved lattice; in the averaged files, the seed mean rounded to a whole number |
| `seed`, `seeds` | The run's random seed; the number of seeds in a mean |
| `step` | Time step |
| `old`, `new` | Seed-mean density, in particles per site, of the cluster (`old`) and of the busiest site's island (`new`). `old` is 0 when no site reaches the cluster cutoff. |
| `island`, `steps` | Which late value (`old`, `new`, `densest`) and the steps it averages: first and last sample step, or the saved lattice used |
| `mean`, `seed_min`, `seed_max` | Late density: mean over seeds, and the smallest and largest seed value |
| `islands` | Number of islands in the lattice |
| `*_mass`, `*_sites` | Particles and sites of the busiest site's island (`new`), the heaviest, the densest, and `densest2`: the densest island of at least two sites, blank when there is none |
