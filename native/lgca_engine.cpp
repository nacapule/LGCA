// native/lgca_engine.cpp — implementation of class Lgca (see lgca_engine.hpp).
//
// Every function below is the reference function of the same name in
// lgca/lgca_clean-1.cpp, with lattice[x][y][k] replaced by occ[(y*W + x)*6 + k]. The
// loops, the order of the random draws and the order of the floating-point operations
// are the reference's, line for line; comments point out the places where that order
// matters. Do not "simplify" an expression here without reading docs/FIDELITY.md.
//
// The exceptions are the research options (fermion_field(), channel_cdf_options()),
// which the reference does not have: they follow the Lab's JavaScript engine
// (collideFermion and collideBoson in the <script id="engineSrc"> block of
// lgca-viz/lgca-lab.html) operation for operation. With the default options the engine
// never enters them.

#include "lgca_engine.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lgca {

namespace {

template <class T> inline T sqr(T a) { return a * a; }

// The six neighbours of a site, in the order of the reference's neighbour_field().
// The order fixes the order of the floating-point sums of the field.
constexpr int NEIGHBOUR[NODES][2] = {{+1, 0}, {0, +1}, {-1, 0},
                                     {0, -1}, {+1, -1}, {-1, +1}};

// Streaming: a particle in channel k at (i, j) moves to (i + dx, j + dy), with the
// offsets of the reference's streaming block (periodic boundaries).
constexpr int STREAM[NODES][2] = {{+1, 0}, {+1, -1}, {0, -1},
                                  {-1, 0}, {-1, +1}, {0, +1}};

}  // namespace

// ----------------------------------------------------------------------------------
//  Construction and initial state
// ----------------------------------------------------------------------------------

const char* params_error(const Params& params) {
    if (params.W < 4 || params.H < 4)
        return "lattice width and height must be at least 4";
    if (!dens_supported(params.model, params.dens))
        return "density outside the range the reference program supports";
    if (!(params.alpha >= 0.0 && params.alpha <= 2.0))   // also refuses NaN
        return "alpha must be a number from 0 to 2";
    return nullptr;
}

Lgca::Lgca(const Params& params)
    : W_(params.W), H_(params.H), model_(params.model), dens_(params.dens),
      kernel_(params.kernel), alpha_(params.alpha), boson_field_(params.boson_field),
      boson_align_(params.boson_align) {
    if (const char* error = params_error(params))
        throw std::invalid_argument(error);

    decide_collision();

    // Channel vectors and kernels exactly as the reference constructor builds them:
    // h = sqrt(3)/2 computed (not a decimal literal), dot products in the same order.
    const double h = std::sqrt(3.0) / 2.0;   // sin(60 deg)
    const double v[NODES][2] = {{ 1.0,  0.0}, { 0.5,  h}, {-0.5,  h},
                                {-1.0,  0.0}, {-0.5, -h}, { 0.5, -h}};
    for (int i = 0; i < NODES; ++i) { c_[i][0] = v[i][0]; c_[i][1] = v[i][1]; }
    for (int i = 0; i < NODES; ++i)
        for (int j = 0; j < NODES; ++j) {
            const double dot = c_[i][0] * c_[j][0] + c_[i][1] * c_[j][1];
            J_[i][j]  = dot;
            J2_[i][j] = dot * dot;
        }

    const std::size_t cells = static_cast<std::size_t>(W_) * H_ * NODES;
    occ_.assign(cells, 0);
    src_.assign(cells, 0);

    rng_.seed_like_reference(params.seed);
    init_lattice();
}

// The Lab tests `kernel == "avg" || (kernel == "power" && alpha === 1)` first, then
// `kernel == "power" && alpha !== 0`, at every collision; kernel and alpha change only
// through the constructor and set_options(), so the outcome is decided there, once.
void Lgca::decide_collision() {
    if (kernel_ == Kernel::Avg || (kernel_ == Kernel::Power && alpha_ == 1.0))
        norm_ = Norm::ByCount;
    else if (kernel_ == Kernel::Power && alpha_ != 0.0)
        norm_ = Norm::ByPower;
    else
        norm_ = Norm::None;
    reference_boson_ = boson_field_ == BosonField::Site &&
                       boson_align_ == BosonAlign::Nematic && norm_ == Norm::None;
}

void Lgca::set_options(Kernel kernel, double alpha, BosonField field, BosonAlign align) {
    if (!(alpha >= 0.0 && alpha <= 2.0))   // also refuses NaN, as params_error() does
        throw std::invalid_argument("alpha must be a number from 0 to 2");
    kernel_ = kernel;
    alpha_ = alpha;
    boson_field_ = field;
    boson_align_ = align;
    decide_collision();
}

// The start of every realization in the reference main():
//   - every channel gets int(dens) particles (boson; the integer part, rounded toward
//     zero, so floor(dens) for dens >= 0); the fermion lattice starts empty whatever
//     dens is;
//   - then one random draw per (i, j, k), i outer, j inner, k innermost, ALWAYS, even
//     when the fill probability frac = dens - int(dens) is 0 or negative: a draw below
//     frac adds one particle.
// N (nuparts) counts int(dens) particles per channel for BOTH models. For the fermion
// model with dens >= 1 those particles are never placed, so N is larger than the
// number of particles on the lattice. That is the reference's behaviour (a quirk, kept;
// docs/FIDELITY.md).
void Lgca::init_lattice() {
    std::int64_t nuparts = std::int64_t(int(dens_)) * NODES * W_ * H_;
    const int32_t base = (model_ == Model::Boson) ? int(dens_) : 0;
    std::fill(occ_.begin(), occ_.end(), base);
    std::fill(src_.begin(), src_.end(), 0);

    const double frac = dens_ - int(dens_);
    for (int i = 0; i < W_; ++i)
        for (int j = 0; j < H_; ++j)
            for (int k = 0; k < NODES; ++k)
                if (rng_.next() < frac) {
                    int32_t& cell = occ_[site(i, j) + k];
                    if (model_ == Model::Boson) cell += 1;   // add_part(vector<int>)
                    else cell = 1;                           // add_part(bitset): set the bit
                    ++nuparts;
                }

    nuparts_ = nuparts;
    t_ = 0;
}

// ----------------------------------------------------------------------------------
//  Time step
// ----------------------------------------------------------------------------------

void Lgca::step(double sens) {
    // (1) Collision, sites in the reference order (i outer, j inner). An empty site
    // draws no random numbers and leaves zeros in src.
    std::fill(src_.begin(), src_.end(), 0);
    for (int i = 0; i < W_; ++i)
        for (int j = 0; j < H_; ++j) {
            const std::size_t s = site(i, j);
            if (is_empty(&occ_[s])) continue;
            if (model_ == Model::Fermion) collide_fermion(i, j, sens, &src_[s]);
            else collide_boson(i, j, sens, &src_[s]);
        }

    // (2) Streaming. Every destination channel is written exactly once, so plain
    // assignment both moves the particles and clears the old lattice.
    for (int i = 0; i < W_; ++i)
        for (int j = 0; j < H_; ++j) {
            const std::size_t from = site(i, j);
            for (int k = 0; k < NODES; ++k) {
                const int x = (i + STREAM[k][0] + W_) % W_;
                const int y = (j + STREAM[k][1] + H_) % H_;
                occ_[site(x, y) + k] = src_[from + k];
            }
        }
    ++t_;
}

// ----------------------------------------------------------------------------------
//  Fermion collision: Metropolis over random channel permutations
// ----------------------------------------------------------------------------------

// Start from the current occupancy; METROPOLIS_STEPS times, propose a uniformly random
// permutation of the channels and accept it if log(u) < E_trial - E. Each proposal
// draws 5 numbers (Fisher-Yates) and the acceptance test 1 more: 72 draws per occupied
// site, whatever is accepted.
void Lgca::collide_fermion(int i, int j, double sens, int32_t* out) {
    double h[NODES];
    collision_field(i, j, h);

    int current[NODES];
    const int32_t* s = &occ_[site(i, j)];
    for (int k = 0; k < NODES; ++k) current[k] = (s[k] != 0) ? 1 : 0;
    double E = align_energy(current, h, sens);

    int trial[NODES];
    for (int step = 0; step < METROPOLIS_STEPS; ++step) {
        random_permutation(trial, current);
        const double E_trial = align_energy(trial, h, sens);
        if (std::log(rng_.next()) < E_trial - E) {
            for (int k = 0; k < NODES; ++k) current[k] = trial[k];
            E = E_trial;
        }
    }
    for (int k = 0; k < NODES; ++k) out[k] = current[k];
}

// The field of a fermion collision at (i, j): the reference's neighbour_field(), or
// fermion_field() when the kernel option divides it.
void Lgca::collision_field(int i, int j, double h[NODES]) const {
    if (norm_ == Norm::None) neighbour_field(i, j, h);   // the reference
    else fermion_field(i, j, h);                           // kernel Avg or Power
}

// E = sens * (sum of h over the occupied channels). The sum skips empty channels,
// as the reference does, rather than adding 0 * h.
double Lgca::align_energy(const int conf[NODES], const double h[NODES], double sens) const {
    double s = 0.0;
    for (int k = 0; k < NODES; ++k)
        if (conf[k]) s += h[k];
    return sens * s;
}

// m_k = the number of particles in channel k summed over the 6 neighbours of (i, j),
// accumulated neighbour by neighbour in NEIGHBOUR order (the Lab's neighbourField()).
void Lgca::neighbour_counts(int i, int j, double m[NODES]) const {
    for (int k = 0; k < NODES; ++k) m[k] = 0.0;
    for (const auto& d : NEIGHBOUR) {
        const int x = (i + d[0] + W_) % W_;
        const int y = (j + d[1] + H_) % H_;   // the reference wraps y with ny (correct)
        const int32_t* n = &occ_[site(x, y)];
        for (int k = 0; k < NODES; ++k) m[k] += n[k];
    }
}

// h_a = sum_k J[a][k] * m_k, with m from neighbour_counts(): the reference's field.
void Lgca::neighbour_field(int i, int j, double h[NODES]) const {
    double m[NODES];
    neighbour_counts(i, j, m);
    for (int a = 0; a < NODES; ++a) {
        double s = 0.0;
        for (int k = 0; k < NODES; ++k) s += J_[a][k] * m[k];
        h[a] = s;
    }
}

// The fermion field with the kernel option (the Lab's collideFermion): the reference's
// field h, then divided by M = sum_k m_k, the number of neighbouring particles (kernel
// Avg, or Power with alpha 1), or by pow(M, alpha) (Power), unless M is 0. As in the
// Lab, M is summed again from m and each h_a is divided, not multiplied by 1/M.
void Lgca::fermion_field(int i, int j, double h[NODES]) const {
    double m[NODES];
    neighbour_counts(i, j, m);
    for (int a = 0; a < NODES; ++a) {
        double s = 0.0;
        for (int k = 0; k < NODES; ++k) s += J_[a][k] * m[k];
        h[a] = s;
    }
    if (norm_ == Norm::ByCount) {
        double M = 0;
        for (int k = 0; k < NODES; ++k) M += m[k];
        if (M != 0)
            for (int a = 0; a < NODES; ++a) h[a] /= M;
    } else if (norm_ == Norm::ByPower) {
        double M = 0;
        for (int k = 0; k < NODES; ++k) M += m[k];
        if (M != 0) {
            const double divisor = std::pow(M, alpha_);
            for (int a = 0; a < NODES; ++a) h[a] /= divisor;
        }
    }
}

// Fisher-Yates shuffle of the channel indices (draws for positions 5, 4, 3, 2, 1),
// then out[perm[k]] = in[k].
void Lgca::random_permutation(int out[NODES], const int in[NODES]) {
    int perm[NODES];
    for (int k = 0; k < NODES; ++k) perm[k] = k;
    for (int k = NODES - 1; k > 0; --k) std::swap(perm[k], perm[rand_int(k + 1)]);
    for (int k = 0; k < NODES; ++k) out[k] = 0;
    for (int k = 0; k < NODES; ++k)
        if (in[k]) out[perm[k]] = 1;
}

// ----------------------------------------------------------------------------------
//  Boson collision: redraw every particle from a Boltzmann distribution
// ----------------------------------------------------------------------------------

// The n particles of the site are removed and each is put back into a channel drawn
// from the site's channel distribution: n draws per site.
void Lgca::collide_boson(int i, int j, double sens, int32_t* out) {
    const int32_t* s = &occ_[site(i, j)];
    int n = 0;
    for (int k = 0; k < NODES; ++k) n += s[k];

    double prob[NODES];
    collision_cdf(i, j, n, sens, prob);
    for (int k = 0; k < NODES; ++k) out[k] = 0;
    for (int p = 0; p < n; ++p) out[sample_channel(prob)] += 1;
}

// The channel CDF of a boson collision at (i, j), whose n particles are about to be
// redrawn: the reference's channel_cdf() with the default options, otherwise
// channel_cdf_options().
void Lgca::collision_cdf(int i, int j, int n, double sens, double prob[NODES]) const {
    if (reference_boson_) channel_cdf(&occ_[site(i, j)], sens, prob);
    else channel_cdf_options(i, j, n, sens, prob);
}

// Cumulative distribution over channels, P(c) proportional to
// exp(sens * sum_k n_k J2[k][c]): the NEMATIC kernel and the site's OWN occupancy only
// (the reference's active boson variant). Operation order as in the reference:
// weights, subtract the maximum, exp and running total, divide, cumulate.
void Lgca::channel_cdf(const int32_t* s, double sens, double prob[NODES]) const {
    for (int c = 0; c < NODES; ++c) {
        double sum = 0.0;
        for (int k = 0; k < NODES; ++k) sum += s[k] * J2_[k][c];
        prob[c] = sens * sum;
    }
    const double mx = *std::max_element(prob, prob + NODES);
    double z = 0.0;
    for (int c = 0; c < NODES; ++c) { prob[c] = std::exp(prob[c] - mx); z += prob[c]; }
    for (int c = 0; c < NODES; ++c) prob[c] /= z;
    for (int c = 1; c < NODES; ++c) prob[c] += prob[c - 1];
}

// The channel CDF with the research options, operation for operation the Lab's
// collideBoson (n = the site's particle count):
//   1. channel counts m and their total numparts: the site's own (Site), the 6
//      neighbours' (Neigh), or the neighbours' with the site's added after them (Both);
//   2. weights w_c = sens * sum_k m_k JT[k][c], JT = J2 (nematic) or J (polar);
//   3. subtract the largest weight (first maximum, as the Lab's loop);
//   4. divide by numparts (ByCount) or by pow(numparts, alpha) (ByPower), unless
//      numparts is 0; exp; running total z;
//   5. divide by z and cumulate.
// For an isolated site under Neigh, numparts = 0 and all weights are 0: every channel
// gets probability 1/6.
void Lgca::channel_cdf_options(int i, int j, int n, double sens, double prob[NODES]) const {
    const int32_t* s = &occ_[site(i, j)];
    const double (*JT)[NODES] = (boson_align_ == BosonAlign::Nematic) ? J2_ : J_;

    double numparts;
    if (boson_field_ == BosonField::Site) {
        numparts = n;
        for (int c = 0; c < NODES; ++c) {
            double sum = 0.0;
            for (int k = 0; k < NODES; ++k) sum += s[k] * JT[k][c];
            prob[c] = sens * sum;
        }
    } else {
        double m[NODES];
        neighbour_counts(i, j, m);
        if (boson_field_ == BosonField::Both)
            for (int k = 0; k < NODES; ++k) m[k] += s[k];
        numparts = 0;
        for (int k = 0; k < NODES; ++k) numparts += m[k];
        for (int c = 0; c < NODES; ++c) {
            double sum = 0.0;
            for (int k = 0; k < NODES; ++k) sum += m[k] * JT[k][c];
            prob[c] = sens * sum;
        }
    }

    double mx = prob[0];
    for (int c = 1; c < NODES; ++c)
        if (prob[c] > mx) mx = prob[c];

    double z = 0.0;
    if (norm_ == Norm::ByCount) {
        for (int c = 0; c < NODES; ++c) {
            double v = prob[c] - mx;
            if (numparts != 0) v /= numparts;
            prob[c] = std::exp(v);
            z += prob[c];
        }
    } else if (norm_ == Norm::ByPower) {
        const double divisor = (numparts == 0) ? 1.0 : std::pow(numparts, alpha_);
        for (int c = 0; c < NODES; ++c) {
            prob[c] = std::exp((prob[c] - mx) / divisor);
            z += prob[c];
        }
    } else {
        for (int c = 0; c < NODES; ++c) {
            prob[c] = std::exp(prob[c] - mx);
            z += prob[c];
        }
    }
    for (int c = 0; c < NODES; ++c) prob[c] /= z;
    for (int c = 1; c < NODES; ++c) prob[c] += prob[c - 1];
}

// One draw r; the first channel with r <= CDF; channel 5 if rounding left the CDF
// just below r.
int Lgca::sample_channel(const double prob[NODES]) {
    const double r = rng_.next();
    for (int c = 0; c < NODES; ++c)
        if (r <= prob[c]) return c;
    return NODES - 1;
}

// ----------------------------------------------------------------------------------
//  Observables (no random draws)
// ----------------------------------------------------------------------------------

// The reference's measure_order(): sums over occupied sites (i outer, j inner, all six
// channels, empty ones included), then the magnitudes divided by N = nuparts.
// For dens = 0, N = 0 and the results are NaN (0/0), as in the reference.
Order Lgca::measure() const {
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0, entropy = 0;
    std::int64_t occupied = 0;
    for (int i = 0; i < W_; ++i)
        for (int j = 0; j < H_; ++j) {
            const int32_t* s = &occ_[site(i, j)];
            if (is_empty(s)) continue;
            ++occupied;
            double nloc = 0;
            for (int k = 0; k < NODES; ++k) {
                const double n = s[k];
                x1 += n * c_[k][0];
                y1 += n * c_[k][1];
                x2 += n * (2 * sqr(c_[k][0]) - 1);     // cos(2 theta_k)
                y2 += n * (2 * c_[k][0] * c_[k][1]);   // sin(2 theta_k)
                nloc += n;
            }
            const double p = nloc / double(nuparts_);
            entropy += p * std::log(p);
        }

    Order out;
    out.polar    = std::sqrt(sqr(x1) + sqr(y1)) / double(nuparts_);
    out.nematic  = std::sqrt(sqr(x2) + sqr(y2)) / double(nuparts_);
    out.entropy  = -entropy;
    out.inv_occ  = occupied > 0 ? 1.0 / double(occupied) : 0.0;
    out.occupied = occupied;
    return out;
}

// The reference's measure_band(): local nematic order of the particles on the four
// rings around each occupied site. One deliberate change (docs/FIDELITY.md, "Deliberate
// deviation"): the reference's channel loop starts at m = 1, skipping channel 0; here it
// counts all six (m = 0). The site weight nsite is summed once per ring, 4 times, and
// divided by 4 below, as in the reference.
double Lgca::band() const {
    double preband = 0;
    for (int i = 0; i < W_; ++i)
        for (int j = 0; j < H_; ++j) {
            if (is_empty(&occ_[site(i, j)])) continue;
            double x3 = 0, y3 = 0, nband = 0, nsite = 0;
            for (int k = 1; k < 5; ++k) {              // rings at distance k
                for (int m = 0; m < NODES; ++m) {      // channels: all six (reference: from 1)
                    const int idx[NODES][2] = {{(i + k) % W_, j},
                                               {(i + k) % W_, (j - k + H_) % H_},
                                               {i, (j - k + H_) % H_},
                                               {(i - k + W_) % W_, j},
                                               {(i - k + W_) % W_, (j + k) % H_},
                                               {i, (j + k) % H_}};
                    for (int d = 0; d < NODES; ++d) {
                        const int occ = occ_[site(idx[d][0], idx[d][1]) + m];
                        x3    += occ * (2 * sqr(c_[d][0]) - 1);
                        y3    += occ * (2 * c_[d][0] * c_[d][1]);
                        nband += occ;
                    }
                    nsite += occ_[site(i, j) + m];
                }
            }
            if (nband > 0) preband += nsite * std::sqrt(sqr(x3) + sqr(y3)) / (4 * double(nband));
        }
    return preband / double(nuparts_);
}

}  // namespace lgca
