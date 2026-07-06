// ============================================================================
//  Lattice-gas cellular automaton (LGCA) for collective cell migration
//  on a 2-D hexagonal lattice (6 velocity channels), with periodic boundaries.
//
//  Each cell carries a velocity equal to one of the 6 hexagonal channel
//  directions. One time step is:
//
//    1. COLLISION  - velocities are stochastically reassigned to favour local
//                    alignment. The strength is the "sensitivity" `sens`
//                    (sens > 0 -> alignment; larger -> stronger ordering).
//                    Particle number is conserved by the collision.
//    2. STREAMING  - every particle hops one site along its velocity channel.
//
//  Two particle statistics, chosen at COMPILE TIME via the MODEL macro below:
//
//    FERMION : at most ONE particle per channel (an exclusion principle, so at
//              most NODES particles per site). State = std::bitset<NODES>.
//              Collision = Metropolis sampling over random channel permutations,
//              driven by a POLAR alignment energy with the 6 neighbours.
//
//    BOSON   : unbounded occupancy per channel (no exclusion). State =
//              std::vector<int>. Collision = redraw every particle's channel
//              from a Boltzmann distribution built from a NEMATIC alignment
//              field (see note in channel_cdf() about locality).
//
//  Build:  g++ -O2 -std=c++17 lgca.cpp -o lgca
//          add  -DCIMG  (and put CImg.h where the include expects it) for the
//          live image output; add  -DHISTOG  to dump a per-site occupancy
//          histogram on the final step.
//
//  Run:    ./lgca  <tsteps> <iters> <sens> <dens> [seed]
//            tsteps : number of time steps per realization
//            iters  : number of independent realizations to average over
//            sens   : alignment sensitivity
//            dens   : initial mean occupancy per channel
//                     FERMION: must satisfy 0 <= dens < 1 (per-channel fill
//                              probability; dens = 1 cannot be represented and
//                              yields an EMPTY lattice -- see initialize()).
//                     BOSON  : any dens >= 0; the integer part fills every
//                              channel, the fractional part is added randomly.
//            seed   : optional unsigned int for reproducible runs
//                     (default: derived from the clock).
//
//  Output (one line): sens*dens  <polar> <err>  <nematic> <err>
//                     <spatial_order> <err>  <band> <err>  <inv_occupancy> <err>
// ============================================================================

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#ifdef CIMG
#include "../CImg/CImg.h"          // live graphics output (optional)
#endif

#include "../rng/WELL1024a.h"      // WELLRNG1024a() -> double in [0,1)

// ----------------------------------------------------------------------------
//  Compile-time configuration
// ----------------------------------------------------------------------------
#define FERMION 0
#define BOSON   1
#define MODEL   BOSON               // <-- set to FERMION or BOSON

constexpr int NODES = 6;            // hexagonal lattice: 6 velocity channels
constexpr int XDIM  = 120;           // lattice width
constexpr int YDIM  = 120;           // lattice height
constexpr int METROPOLIS_STEPS = 12;// collision relaxation sweeps (fermion)

template <class T> inline T sqr(T a) { return a * a; }

// Uniform integer in [0, n) from the WELL generator.
inline int rand_int(int n) { return static_cast<int>(WELLRNG1024a() * n); }

// ============================================================================
//  LGCA lattice and collision operators
// ============================================================================
template <class State>
class LGCA {
public:
    int nx, ny;
    std::vector<std::vector<State>> lattice;        // pre-collision configuration
    std::vector<std::vector<State>> lattice_temp;   // post-collision, pre-streaming
    double c[NODES][2];                             // channel velocity unit vectors
    double J[NODES][NODES];                         // c_i . c_j        (polar kernel)
    double J2[NODES][NODES];                        // (c_i . c_j)^2    (nematic kernel)
    std::vector<double> prob;                       // scratch: channel CDF (boson)

    LGCA(int nx_, int ny_) : nx(nx_), ny(ny_), prob(NODES, 0.0) {
        const double h = std::sqrt(3.0) / 2.0;      // sin(60 deg)
        const double v[NODES][2] = {{ 1.0,  0.0}, { 0.5,  h}, {-0.5,  h},
                                    {-1.0,  0.0}, {-0.5, -h}, { 0.5, -h}};
        for (int i = 0; i < NODES; ++i) { c[i][0] = v[i][0]; c[i][1] = v[i][1]; }

        for (int i = 0; i < NODES; ++i)
            for (int j = 0; j < NODES; ++j) {
                double dot = c[i][0] * c[j][0] + c[i][1] * c[j][1];
                J[i][j]  = dot;
                J2[i][j] = dot * dot;
            }

        lattice.assign(nx, std::vector<State>(ny));
        lattice_temp.assign(nx, std::vector<State>(ny));
    }

    // ----- FERMION collision -------------------------------------------------
    // Metropolis sampling: start from the current occupancy, repeatedly propose
    // a random permutation of the channels, and accept it with the usual
    // Metropolis rule for the polar alignment energy E = sens * sum_i n_i h_i.
    void collide(std::bitset<NODES>& out, int x, int y, double sens) {
        double h[NODES];
        neighbour_field(x, y, h);                   // h_i from the 6 neighbours
        out = lattice[x][y];
        double E = align_energy(out, h, sens);

        std::bitset<NODES> trial;
        for (int k = 0; k < METROPOLIS_STEPS; ++k) {
            random_permutation(trial, out);
            double E_trial = align_energy(trial, h, sens);
            // accept with probability min(1, exp(E_trial - E))
            if (std::log(WELLRNG1024a()) < E_trial - E) {
                out = trial;
                E   = E_trial;
            }
        }
    }

    // ----- BOSON collision ---------------------------------------------------
    // Redraw each particle's channel independently from the Boltzmann
    // distribution computed in channel_cdf().
    void collide(std::vector<int>& out, int x, int y, double sens) {
        int n = 0;
        for (int i = 0; i < NODES; ++i) n += lattice[x][y][i];
        channel_cdf(x, y, sens);
        out.assign(NODES, 0);
        for (int p = 0; p < n; ++p) out[sample_channel()] += 1;
    }

private:
    // Polar alignment energy of a (fermion) configuration in field h.
    double align_energy(const std::bitset<NODES>& conf,
                        const double h[NODES], double sens) const {
        double s = 0.0;
        for (int i = 0; i < NODES; ++i)
            if (conf[i]) s += h[i];
        return sens * s;
    }

    // h_i = sum_j J[i][j] * m_j, where m_j is the total occupancy of channel j
    // summed over the 6 hexagonal neighbours. This is algebraically identical
    // to the original triple loop but evaluates the energy in O(NODES) per
    // Metropolis proposal instead of O(NODES^2 * neighbours).
    void neighbour_field(int x, int y, double h[NODES]) const {
        static const int nb[6][2] = {{+1, 0}, {0, +1}, {-1, 0},
                                     {0, -1}, {+1, -1}, {-1, +1}};
        double m[NODES] = {0};
        for (const auto& d : nb) {
            int X = (x + d[0] + nx) % nx;
            int Y = (y + d[1] + ny) % ny;           // NOTE: uses ny (see report)
            for (int j = 0; j < NODES; ++j) m[j] += lattice[X][Y][j];
        }
        for (int i = 0; i < NODES; ++i) {
            double s = 0.0;
            for (int j = 0; j < NODES; ++j) s += J[i][j] * m[j];
            h[i] = s;
        }
    }

    // Uniformly random permutation of channel occupancies (Fisher-Yates).
    // Conserves particle number and is a symmetric Metropolis proposal.
    void random_permutation(std::bitset<NODES>& out,
                            const std::bitset<NODES>& in) const {
        int perm[NODES];
        for (int i = 0; i < NODES; ++i) perm[i] = i;
        for (int i = NODES - 1; i > 0; --i) std::swap(perm[i], perm[rand_int(i + 1)]);
        out.reset();
        for (int i = 0; i < NODES; ++i)
            if (in[i]) out.set(perm[i]);
    }

    // Build the per-channel Boltzmann CDF for the boson collision.
    // The alignment field here is NEMATIC and uses ONLY the local site's
    // occupancy (no neighbours). This matches the active version of the
    // original code; the neighbour-coupled variant was commented out there.
    // Verify this locality is what you intend (the fermion model, by contrast,
    // aligns with neighbours).
    void channel_cdf(int x, int y, double sens) {
        const State& site = lattice[x][y];
        for (int j = 0; j < NODES; ++j) {
            double s = 0.0;
            for (int i = 0; i < NODES; ++i) s += site[i] * J2[i][j];
            prob[j] = sens * s;
        }
        // subtract the max for numerical stability (does not change the result)
        double mx = *std::max_element(prob.begin(), prob.end());
        double z = 0.0;
        for (int j = 0; j < NODES; ++j) { prob[j] = std::exp(prob[j] - mx); z += prob[j]; }
        for (int j = 0; j < NODES; ++j) prob[j] /= z;
        for (int j = 1; j < NODES; ++j) prob[j] += prob[j - 1];   // -> CDF
    }

    // Draw a channel index from the CDF in prob[].
    int sample_channel() const {
        double r = WELLRNG1024a();
        for (int i = 0; i < NODES; ++i)
            if (r <= prob[i]) return i;
        return NODES - 1;            // FIX: guard the floating-point CDF edge case
    }
};

// ============================================================================
//  State-type helpers (overloaded for fermion / boson)
// ============================================================================
inline bool is_empty(const std::bitset<NODES>& site) { return site.none(); }
inline bool is_empty(const std::vector<int>& site) {
    for (int v : site) if (v > 0) return false;
    return true;
}

inline int site_count(const std::bitset<NODES>& site) { return static_cast<int>(site.count()); }
inline int site_count(const std::vector<int>& site) {
    int n = 0; for (int v : site) n += v; return n;
}

// FERMION init: each channel occupied independently with probability `fill`
// (so `dens` is a per-channel probability and must lie in [0,1)).
inline void initialize(std::bitset<NODES>& conf, double /*fill ignored here*/ = 0.0) {
    conf.reset();
}
// BOSON init: put floor(fill) particles in every channel.
inline void initialize(std::vector<int>& conf, double fill = 0.0) {
    conf.assign(NODES, static_cast<int>(fill));
}

inline void add_part(std::bitset<NODES>& conf, int chan) { conf.set(chan); }
inline void add_part(std::vector<int>& conf, int chan) { conf[chan] += 1; }

// ============================================================================
//  Observables (computed each step, accumulated across realizations)
// ============================================================================
template <class State>
struct StepStats {
    double polar    = 0;   // |sum of velocity vectors| / N
    double nematic  = 0;   // |sum of nematic (doubled-angle) vectors| / N
    double entropy  = 0;   // Shannon entropy of the spatial occupancy
    double band     = 0;   // local nematic-band order parameter
    double inv_occ  = 0;   // 1 / (number of occupied sites)
};

// Polar + nematic order parameters and spatial entropy.
template <class State>
void measure_order(const LGCA<State>& s, long nuparts, StepStats<State>& out) {
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0, entropy = 0;
    long occupied = 0;
    for (int i = 0; i < s.nx; ++i)
        for (int j = 0; j < s.ny; ++j) {
            const State& site = s.lattice[i][j];
            if (is_empty(site)) continue;
            ++occupied;
            double nloc = 0;
            for (int k = 0; k < NODES; ++k) {
                double n = site[k];
                x1 += n * s.c[k][0];
                y1 += n * s.c[k][1];
                x2 += n * (2 * sqr(s.c[k][0]) - 1);   // cos(2 theta_k)
                y2 += n * (2 * s.c[k][0] * s.c[k][1]); // sin(2 theta_k)
                nloc += n;
            }
            double p = nloc / double(nuparts);
            entropy += p * std::log(p);
        }
    out.polar    = std::sqrt(sqr(x1) + sqr(y1)) / double(nuparts);
    out.nematic  = std::sqrt(sqr(x2) + sqr(y2)) / double(nuparts);
    out.entropy  = -entropy;
    out.inv_occ  = occupied > 0 ? 1.0 / double(occupied) : 0.0;
}

// Local nematic "band" order parameter, averaged over occupied sites.
// Kept arithmetically identical to the original implementation.
// NOTE: this scan starts the channel loop at k = 1, i.e. it SKIPS channel 0,
// unlike measure_order() which includes all channels. Verify this is intended.
template <class State>
double measure_band(const LGCA<State>& s, long nuparts) {
    double preband = 0;
    for (int i = 0; i < s.nx; ++i)
        for (int j = 0; j < s.ny; ++j) {
            if (is_empty(s.lattice[i][j])) continue;
            double x3 = 0, y3 = 0, nband = 0, nsite = 0;
            for (int k = 1; k < 5; ++k) {                 // rings at distance k
                for (int m = 1; m < NODES; ++m) {         // channels (skips 0)
                    const int idx[6][2] = {{(i + k) % s.nx, j},
                                           {(i + k) % s.nx, (j - k + s.ny) % s.ny},
                                           {i, (j - k + s.ny) % s.ny},
                                           {(i - k + s.nx) % s.nx, j},
                                           {(i - k + s.nx) % s.nx, (j + k) % s.ny},
                                           {i, (j + k) % s.ny}};
                    for (int d = 0; d < NODES; ++d) {
                        int occ = s.lattice[idx[d][0]][idx[d][1]][m];
                        x3    += occ * (2 * sqr(s.c[d][0]) - 1);
                        y3    += occ * (2 * s.c[d][0] * s.c[d][1]);
                        nband += occ;
                    }
                    nsite += s.lattice[i][j][m];
                }
            }
            if (nband > 0) preband += nsite * std::sqrt(sqr(x3) + sqr(y3)) / (4 * double(nband));
        }
    return preband / double(nuparts);
}

// ============================================================================
//  main
// ============================================================================
int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "usage: " << argv[0]
                  << " <tsteps> <iters> <sens> <dens> [seed]\n";
        return 1;
    }
    const int    tsteps = std::atoi(argv[1]);
    const int    iters  = std::atoi(argv[2]);
    const double sens   = std::atof(argv[3]);
    const double dens   = std::atof(argv[4]);

    // --- seed the WELL1024a generator ---------------------------------------
    // Seed all 32 state words with DISTINCT values (a one-line LCG mixes them),
    // so two runs launched in the same second do not share a stream. Pass an
    // explicit seed as argv[5] for reproducible runs.
    unsigned int seed = (argc > 5) ? static_cast<unsigned int>(std::strtoul(argv[5], nullptr, 10))
                                   : static_cast<unsigned int>(std::time(nullptr));
    unsigned int init[32];
    for (int i = 0; i < 32; ++i) { init[i] = seed; seed = seed * 1664525u + 1013904223u; }
    InitWELLRNG1024a(init);

#if MODEL == FERMION
    using State = std::bitset<NODES>;
#else
    using State = std::vector<int>;
#endif

    LGCA<State> sim(XDIM, YDIM);
    State temp;

#ifdef CIMG
    cimg_library::CImg<unsigned char> frame(4 * XDIM + 4, 4 * YDIM + 4, 1, 3, 0);
    cimg_library::CImgDisplay display(frame, "result");
    frame.fill(0);
    unsigned char rgb[3] = {255, 255, 255};
#endif

    // --- time-series accumulators (sum and sum of squares over realizations) -
    std::vector<double> polar(tsteps, 0),   polar2(tsteps, 0);
    std::vector<double> nemat(tsteps, 0),   nemat2(tsteps, 0);
    std::vector<double> entro(tsteps, 0),   entro2(tsteps, 0);
    std::vector<double> band(tsteps, 0),    band2(tsteps, 0);
    std::vector<double> invocc(tsteps, 0),  invocc2(tsteps, 0);
    std::vector<double> polar_max(tsteps, 0), nemat_max(tsteps, 0);
    std::vector<double> entro_min(tsteps, 0), band_max(tsteps, 0);

    // ========================== realizations ================================
    for (int r = 0; r < iters; ++r) {

        // particle bookkeeping: integer part fills every channel uniformly
        long nuparts = long(int(dens)) * NODES * XDIM * YDIM;

        // initialize the lattice
        for (int i = 0; i < XDIM; ++i)
            for (int j = 0; j < YDIM; ++j) {
                initialize(sim.lattice[i][j], dens);   // floor(dens) per channel (boson)
                initialize(sim.lattice_temp[i][j]);
            }
        // add the fractional part: each channel gets one more particle with
        // probability frac(dens). For the fermion model this is the ONLY source
        // of particles (so dens is the per-channel occupation probability).
        double frac = dens - int(dens);
        for (int i = 0; i < XDIM; ++i)
            for (int j = 0; j < YDIM; ++j)
                for (int k = 0; k < NODES; ++k)
                    if (WELLRNG1024a() < frac) { add_part(sim.lattice[i][j], k); ++nuparts; }

        // ============================ time loop =============================
        for (int t = 0; t < tsteps; ++t) {

            // (1) COLLISION: choose a post-collision configuration per site.
            for (int i = 0; i < XDIM; ++i)
                for (int j = 0; j < YDIM; ++j)
                    if (!is_empty(sim.lattice[i][j])) {
                        sim.collide(temp, i, j, sens);
                        sim.lattice_temp[i][j] = temp;
                    }

            // (2) STREAMING: each channel hops to its neighbour, then clear temp.
            // Every destination (cell, channel) is written exactly once, so the
            // plain assignment both moves particles and erases the old lattice.
            for (int i = 0; i < XDIM; ++i)
                for (int j = 0; j < YDIM; ++j) {
                    sim.lattice[(i + 1) % XDIM][j][0]                       = sim.lattice_temp[i][j][0];
                    sim.lattice[(i + 1) % XDIM][(j - 1 + YDIM) % YDIM][1]   = sim.lattice_temp[i][j][1];
                    sim.lattice[i][(j - 1 + YDIM) % YDIM][2]                = sim.lattice_temp[i][j][2];
                    sim.lattice[(i - 1 + XDIM) % XDIM][j][3]                = sim.lattice_temp[i][j][3];
                    sim.lattice[(i - 1 + XDIM) % XDIM][(j + 1) % YDIM][4]   = sim.lattice_temp[i][j][4];
                    sim.lattice[i][(j + 1) % YDIM][5]                       = sim.lattice_temp[i][j][5];
                    for (int k = 0; k < NODES; ++k) sim.lattice_temp[i][j][k] = 0;
                }

#ifdef HISTOG
            if (t == tsteps - 1) {
                std::stringstream fname;
                fname << "histogram_t" << t << "_s" << sens << "_d" << dens << ".csv";
                std::ofstream hist(fname.str(), std::ios::out | std::ios::app);
                for (int i = 0; i < XDIM; ++i)
                    for (int j = 0; j < YDIM; ++j)
                        if (!is_empty(sim.lattice[i][j]))
                            hist << double(site_count(sim.lattice[i][j])) / double(NODES) << ", ";
                hist.close();
            }
#endif

            // (3) OBSERVABLES.
            StepStats<State> stats;
            measure_order(sim, nuparts, stats);
            stats.band = measure_band(sim, nuparts);

            polar[t]  += stats.polar;    polar2[t]  += sqr(stats.polar);
            nemat[t]  += stats.nematic;  nemat2[t]  += sqr(stats.nematic);
            entro[t]  += stats.entropy;  entro2[t]  += sqr(stats.entropy);
            band[t]   += stats.band;     band2[t]   += sqr(stats.band);
            invocc[t] += stats.inv_occ;  invocc2[t] += sqr(stats.inv_occ);

            if (stats.polar   > polar_max[t]) polar_max[t] = stats.polar;
            if (stats.nematic > nemat_max[t]) nemat_max[t] = stats.nematic;
            if (stats.band    > band_max[t])  band_max[t]  = stats.band;
            if (r == 0 || stats.entropy < entro_min[t]) entro_min[t] = stats.entropy;

            // (4) OPTIONAL live image (colour each site by its occupancy).
#ifdef CIMG
            frame.fill(0);
            for (int i = 0; i < XDIM; ++i)
                for (int j = 0; j < YDIM; ++j) {
                    if (is_empty(sim.lattice[i][j])) { rgb[0] = rgb[1] = rgb[2] = 0; }
                    else {
                        double brightness = site_count(sim.lattice[i][j]) / double(NODES * dens) * 100 + 155;
                        rgb[0] = rgb[1] = static_cast<unsigned char>(brightness > 510 ? 255 : brightness / 2.0);
                        rgb[2]          = static_cast<unsigned char>(brightness > 255 ? 255 : brightness);
                    }
                    for (int ip = 0; ip < 4; ++ip)
                        for (int jp = 0; jp < 4; ++jp)
                            frame.draw_point(4 * ((i + (j / 2)) % XDIM) + ip + 2 * (j % 2),
                                             4 * j + jp, rgb);
                }
            display.display(frame);
#endif
        } // time loop
    } // realizations

    // --- output: averages and standard errors at the final time step --------
    auto stderr_of = [iters](double sum, double sum2) {
        return std::sqrt(std::fabs(sqr(sum / iters) - sum2 / iters) / iters);
    };
    int T = tsteps - 1;
    std::cout << sens * dens << " "
              << polar[T] / iters  << " " << stderr_of(polar[T],  polar2[T])  << " "
              << nemat[T] / iters  << " " << stderr_of(nemat[T],  nemat2[T])  << " "
              << 1.0 - (entro[T] / iters) / std::log(double(XDIM * YDIM)) << " "
              << stderr_of(entro[T], entro2[T]) / std::log(double(XDIM * YDIM)) << " "
              << band[T] / iters   << " " << stderr_of(band[T],   band2[T])   << " "
              << invocc[T] / iters << " " << stderr_of(invocc[T], invocc2[T]) << "\n";

    return 0;
}
