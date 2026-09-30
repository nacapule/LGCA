# Model and measurements

The lattice has W × H sites, periodic boundaries and six velocity channels.
Channel k points along c_k = (cos(kπ/3), sin(kπ/3)), for k = 0,…,5.
Each step first changes particle directions at each occupied site, conserving
that site's particle count, then streams each particle to its neighbour.
The total number of particles N is conserved.

Fermions allow at most one particle per channel. Their collision uses the
six-neighbour polar field and 12 Metropolis proposals that permute the channels;
a proposal is accepted when log(u) < E_trial − E, with
E = sensitivity × the sum of alignment fields in occupied channels.
Bosons allow any number of particles per channel. Each particle independently
redraws its direction from the same collision probability at its site.

The C++ boson default uses the site's own population, nematic alignment
J_kc = (c_k · c_c)^2 and an additive field.
The September 30 study uses bosons, polar alignment J_kc = c_k · c_c,
the centre plus its six neighbours, and sensitivity s = 6.

## Alpha kernel

Let m_k be the particles in channel k over the contributing sites and M = Σ_k m_k.
For the study, this includes the centre and all six neighbours, including the
particle being redrawn. Define w_c = s Σ_k m_k J_kc. The collision probability is

    p_c = exp((w_c − max_d w_d) / M^alpha)
          / Σ_d exp((w_d − max_e w_e) / M^alpha).

When M = 0, all six directions have probability 1/6.
Alpha 0 is the additive sum; alpha 1 is the particle-count average.
For fermions, the neighbour field is divided by M^alpha before the Metropolis test.
At fixed directional composition, the field strength scales as M^(1−alpha).

## Observables used

The input `dens` is occupancy per channel; expected initial particles per site
are 6 × dens. The actual average density is rho = N/(W × H).
Boson initialization assigns floor(dens) particles to each channel and adds
one with probability dens − floor(dens).

Let n_x = Σ_k n_xk. Qualifying cluster sites satisfy
n_x ≥ max(1, ceil(2 × rho)). Join qualifying sites through the six-neighbour
graph with periodic boundaries. The largest cluster is the component with
the most particles, not the most sites.

Cluster density = cluster mass / number of cluster sites, in particles/site.
It is zero if no sites qualify.

Condensate size N_c = max_(x,k) n_xk, stored as `kmax`: the largest occupancy
of a single velocity channel anywhere in the box. It is distinct from the
largest cluster's mass and from the maximum total population of one site.
The size comparison increases W = H at fixed dens, so N increases with box area.

Late values are averages over steps 15,100–20,000 (the last quarter of the run),
first per run, then over seeds.
