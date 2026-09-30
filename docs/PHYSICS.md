# Model and what I measure

Hexagonal lattice W × H with periodic borders, 6 velocity channels per site. Channel k
points in direction c_k = (cos(kπ/3), sin(kπ/3)), k = 0,…,5. One step is: first at each
occupied site the particles change direction (the number of particles in the site
stays the same), then every particle moves one site in its direction. Total number of
particles N is conserved.

Fermions: max one particle per channel. The collision uses the polar field of the 6
neighbours and 12 Metropolis proposals that permute the channels; a proposal is
accepted if log(u) < E_trial − E, with E = sensitivity × sum of the alignment field
over the occupied channels.

Bosons: any number of particles per channel. Each particle of the site picks again its
direction, independently, from the same collision probability.

In the C++ the boson default uses only the site's own particles, nematic alignment
J_kc = (c_k · c_c)^2 and the field as a sum. For the Sep 30 study I use bosons with
polar alignment J_kc = c_k · c_c, field from the centre site plus its 6 neighbours, and
sensitivity s = 6.

## Alpha

m_k = particles in channel k in the sites that contribute, M = Σ_k m_k. In the study
this is the centre and the 6 neighbours, counting also the particle that is choosing.
With w_c = s Σ_k m_k J_kc the probability to go to channel c is

    p_c = exp((w_c − max_d w_d) / M^alpha)
          / Σ_d exp((w_d − max_e w_e) / M^alpha).

If M = 0 all 6 directions have probability 1/6. Alpha 0 is the sum, alpha 1 is the
average per particle. For fermions the neighbour field is divided by M^alpha before the
Metropolis test. With fixed proportions of directions, the field grows like M^(1−alpha).

## What I measure

The input `dens` is particles per channel, so at the start there are 6 × dens particles
per site on average. The real average density is rho = N/(W × H). For bosons each
channel starts with floor(dens) particles, plus one more with probability
dens − floor(dens).

Cluster: n_x = Σ_k n_xk is the particles in site x. A site counts if
n_x ≥ max(1, ceil(2 × rho)), so at least 2 times the average. Sites that count and
touch (6 neighbours, periodic) make one cluster. The largest cluster is the one with most
particles, not with most sites.

Cluster density = particles of the cluster / number of sites of the cluster, in
particles per site (0 if no site counts).

Condensate size N_c = max over x, k of n_xk (`kmax` in the data): the most particles in
one channel of one site in all the box, like in the 2023 paper. It is not the same as
the mass of the cluster or the max particles in one site. For the box size runs I make
W = H bigger with the same dens, so N grows with the area.

Late values are the mean over steps 15,100 to 20,000 (last quarter of the run), first
for each run and then over the seeds.
