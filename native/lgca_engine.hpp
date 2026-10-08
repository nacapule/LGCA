// native/lgca_engine.hpp — the hexagonal lattice-gas engine as a C++ class.
//
// A port of the reference program lgca/lgca_clean-1.cpp with the lattice size and the
// model chosen at run time instead of at compile time. With the same seed it produces
// the same lattice states and the same random-number stream as the reference, bit for
// bit (checked against the reference program). The rules that keep it that way are in
// docs/FIDELITY.md:
//   - the same random draws in the same order (initial fill, collisions);
//   - the same floating-point operations in the same order (energies, Boltzmann
//     weights, observables), compiled with -ffp-contract=off and never -ffast-math;
//   - the reference's quirks kept as they are (see measure() below, and init_lattice()
//     for the fermion model with dens >= 1); one deliberate exception, band() counts
//     channel 0 (docs/FIDELITY.md, "The band correction").
// Some loops are written differently from the reference for speed (streaming, the
// fermion occupancy as a bit mask, the boson channel draw, cached powers); each says in
// lgca_engine.cpp why it gives the same values, draws and states.
//
// It also has the Lab's research options (kernel sum / avg / power with exponent alpha,
// boson field from the site / its neighbours / both, nematic or polar boson alignment;
// docs/PHYSICS.md). With their defaults the engine is the reference program; with other
// values it reproduces the Lab's JavaScript engine (lgca-viz/lgca-lab.html, the
// <script id="engineSrc"> block), whose arithmetic the option code follows line by line.
//
// Memory layout, shared with the Lab's JavaScript engine (lgca-viz/lgca-lab.html) so a
// later WebAssembly build can hand the arrays to the page unchanged:
//   occ[(j*W + i)*6 + k] = number of particles at site (i, j) moving in channel k,
//   i = x in [0, W), j = y in [0, H), k = 0..5. Fermion values are 0 or 1.
// The reference stores the same numbers as lattice[i][j][k].

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "well1024a.hpp"

namespace lgca {

constexpr int NODES = 6;              // velocity channels of the hexagonal lattice
constexpr int METROPOLIS_STEPS = 12;  // fermion collision: proposals per site

// The power kernel divides by pow(M, alpha) for a particle count M. The engine keeps the
// values for M = 1 .. POW_TABLE_SIZE - 1 in a table it allocates once, at construction
// (POW_TABLE_SIZE doubles, 512 KiB); larger counts call std::pow every time. A caller that
// checks memory before constructing (native/wasm/lgca_wasm.cpp) counts this table too.
constexpr std::size_t POW_TABLE_SIZE = std::size_t(1) << 16;

enum class Model { Fermion, Boson };

// The Lab's research options (docs/PHYSICS.md). The first value of each is the default
// and the reference program's behaviour. Names as in the Lab's engine: kernel, alpha,
// bosonField, bosonAlign.
//
// Kernel: what the field is divided by before it enters the collision.
//   Sum    nothing (the reference);
//   Avg    the number of particles M that built the field (nothing if M = 0);
//   Power  M^alpha (nothing if M = 0). alpha = 0 runs the Sum arithmetic and alpha = 1
//          the Avg arithmetic, exactly as the Lab does, so those two endpoints give the
//          same trajectories as Sum and Avg.
// Fermion model: the field is always built from the 6 neighbours with the polar kernel;
// BosonField and BosonAlign do not apply to it (the Lab ignores them too).
enum class Kernel { Sum, Avg, Power };
// Boson model: whose particles build a site's field.
//   Site   the site's own particles (the reference);
//   Neigh  the 6 neighbouring sites, not the site itself;
//   Both   the 6 neighbours plus the site (seven sites).
enum class BosonField { Site, Neigh, Both };
// Boson model: the alignment kernel, nematic (c_i . c_j)^2 (the reference) or polar
// c_i . c_j.
enum class BosonAlign { Nematic, Polar };

// The densities the engine accepts when a realization runs.
//   - The reference converts dens to an int, which is defined only if its integer part
//     fits: -2147483649 < dens < 2147483648 (and dens is a number). Outside that range
//     the reference always has undefined behaviour.
//   - Boson model: every channel starts with int(dens) particles, plus one with
//     probability dens - int(dens), and the collision sums a site's 6 channels in an
//     int. Above (2^31 - 1) / 6 = 357913941 that sum CAN overflow in the very first
//     step: whether it does depends on the random fill (a seed that adds at most one
//     particle to every site stays defined). Refusing every boson dens above that value
//     is therefore conservative, but it guarantees the first collision is defined.
// Negative densities are defined: int(dens) rounds toward zero, so the channels start at
// int(dens) <= 0 and, as dens - int(dens) <= 0, nothing is added; the reference counts
// that as an empty lattice.
constexpr double MAX_BOSON_DENS = 357913941.0;
inline bool dens_supported(Model model, double dens) {
    if (!(dens > -2147483649.0 && dens < 2147483648.0)) return false;   // also NaN
    return model == Model::Fermion || dens <= MAX_BOSON_DENS;
}

// Everything the reference fixes at compile time (size, model) or reads once in main()
// (density, seed), and the Lab's research options. The sensitivity is not here: it is
// an argument of step(), as in the reference, where it can change between steps.
struct Params {
    int W = 120;                  // lattice width  (XDIM in the reference)
    int H = 120;                  // lattice height (YDIM in the reference)
    Model model = Model::Boson;   // the reference's default build is BOSON
    double dens = 0.0;            // particles per channel (integer part + fill probability)
    uint32_t seed = 0;            // seeds WELL1024a exactly as the reference main() does

    // Research options; the defaults are the reference program.
    Kernel kernel = Kernel::Sum;
    double alpha = 0.0;           // exponent of the Power kernel, 0 <= alpha <= 2
    BosonField boson_field = BosonField::Site;
    BosonAlign boson_align = BosonAlign::Nematic;
};

// Why the constructor would refuse these parameters (the message it throws), or nullptr
// if it accepts them. For callers that cannot catch exceptions (the WebAssembly build,
// native/wasm/lgca_wasm.cpp, checks first and never constructs a refused engine).
const char* params_error(const Params& params);

// The observables of measure_order() in the reference.
struct Order {
    double polar = 0.0;     // |sum of particle velocities| / N
    double nematic = 0.0;   // |sum of doubled-angle vectors| / N
    double entropy = 0.0;   // Shannon entropy of where the particles sit
    double inv_occ = 0.0;   // 1 / (number of occupied sites), 0 if none
    std::int64_t occupied = 0;   // number of occupied sites
};

class Lgca {
public:
    // Seeds the generator and fills the lattice for the first realization (the same
    // random draws as the first realization of the reference program). Throws
    // std::invalid_argument if W or H is below 4 (the reference's band scan looks 4
    // sites away and can read outside a smaller lattice), if dens_supported() is false
    // (see above), or if alpha is not a number from 0 to 2 (the Lab's range; the Lab
    // turns a non-number into 0, this engine refuses it).
    explicit Lgca(const Params& params);

    // Start a new realization: refill the lattice from the current generator state.
    // The reference does this at the start of every realization without reseeding, so
    // realizations 2, 3, ... continue the same random stream.
    void init_lattice();

    // One time step: collision at every occupied site, then streaming.
    void step(double sens);

    // Observables of the current lattice. They draw no random numbers.
    Order measure() const;
    double band() const;

    // Size, model and bookkeeping.
    int W() const { return W_; }
    int H() const { return H_; }
    Model model() const { return model_; }
    double dens() const { return dens_; }
    // Counters are 64-bit on every platform (`long` is only 32 bits in WebAssembly, where
    // N = int(dens) * 6 * W * H could overflow for a large dens).
    std::int64_t nuparts() const { return nuparts_; }   // N, the particle count the observables divide by
    std::int64_t t() const { return t_; }               // steps since the last init_lattice()
    Kernel kernel() const { return kernel_; }
    double alpha() const { return alpha_; }
    BosonField boson_field() const { return boson_field_; }
    BosonAlign boson_align() const { return boson_align_; }

    // The lattice (occ) and the post-collision lattice of the last step (src), both
    // W*H*6 int32 values in the layout described at the top of this file.
    int32_t* occ() { return occ_.data(); }
    const int32_t* occ() const { return occ_.data(); }
    const int32_t* src() const { return src_.data(); }
    std::size_t cell_count() const { return occ_.size(); }

    // Generator state in the Lab's layout: 32 state words, then the index.
    void get_rng_state(uint32_t out[33]) const { rng_.get_state(out); }
    void set_rng_state(const uint32_t in[33]) { rng_.set_state(in); }

    // For replaying a saved state into an engine (with occ() and set_rng_state()): N and
    // the step counter of that state. N must be the saved value, not a recount of the
    // lattice (for the fermion model with dens >= 1 they differ; see init_lattice()).
    // Neither enters the dynamics; N divides the observables.
    void set_nuparts(std::int64_t n) { nuparts_ = n; }
    void set_t(std::int64_t t) { t_ = t; }

    // Change the research options between steps, as the Lab does when its user switches
    // them during a run (the Lab's engine reads them at every collision, so a change acts
    // from the next step on, and so it does here). Throws std::invalid_argument, and
    // changes nothing, if alpha is not a number from 0 to 2.
    void set_options(Kernel kernel, double alpha, BosonField field, BosonAlign align);

private:
    // A test program compares the private collision arithmetic below (channel CDF,
    // neighbour field, alignment energy, velocity tables) value by value with the
    // reference and, for the research options, with the Lab's engine; it defines this
    // struct to reach them.
    friend struct LgcaProbe;

    // How the kernel option divides the field, decided once in the constructor with the
    // Lab's test, which it repeats at every collision:
    //   ByCount  kernel Avg, or Power with alpha == 1: divide by M (the Avg arithmetic);
    //   ByPower  kernel Power with alpha other than 0 and 1: divide by pow(M, alpha);
    //   None     kernel Sum, or Power with alpha == 0 (the Sum arithmetic).
    enum class Norm { None, ByCount, ByPower };
    // Sets norm_ and reference_boson_ from the four options (constructor, set_options).
    void decide_collision();

    // Offset of site (i, j) in occ / src.
    std::size_t site(int i, int j) const {
        return (static_cast<std::size_t>(j) * W_ + i) * NODES;
    }
    // The reference's is_empty(): no channel holds a particle.
    static bool is_empty(const int32_t* s) {
        for (int k = 0; k < NODES; ++k)
            if (s[k] > 0) return false;
        return true;
    }
    // The reference's rand_int(n): uniform integer in [0, n).
    int rand_int(int n) { return static_cast<int>(rng_.next() * n); }

    void collide_fermion(int i, int j, double sens, int32_t* out);
    void collide_boson(int i, int j, double sens, int32_t* out);
    void collision_field(int i, int j, double h[NODES]) const;
    void collision_cdf(int i, int j, int n, double sens, double prob[NODES]) const;
    void neighbour_counts(int i, int j, double m[NODES]) const;
    void neighbour_field(int i, int j, double h[NODES]) const;
    void fermion_field(int i, int j, double h[NODES]) const;
    // A fermion site's occupancy as a bit mask: bit k set = channel k holds a particle.
    double align_energy(unsigned conf, const double h[NODES], double sens) const;
    double align_energy(const int conf[NODES], const double h[NODES], double sens) const;
    void channel_cdf(const int32_t* s, double sens, double prob[NODES]) const;
    void channel_cdf_options(int i, int j, int n, double sens, double prob[NODES]) const;
    int sample_channel(const double prob[NODES]);
    double pow_alpha(double m) const;

    int W_, H_;
    Model model_;
    double dens_;
    std::int64_t nuparts_ = 0;
    std::int64_t t_ = 0;

    Kernel kernel_;
    double alpha_;
    BosonField boson_field_;
    BosonAlign boson_align_;
    Norm norm_;                 // the kernel's division, from kernel_ and alpha_
    bool reference_boson_;      // boson collision is the reference's (site, nematic, no division)

    double c_[NODES][2];        // channel velocity unit vectors
    double J_[NODES][NODES];    // c_i . c_j       (polar kernel, fermion field)
    double J2_[NODES][NODES];   // (c_i . c_j)^2   (nematic kernel, boson weights)

    std::vector<int32_t> occ_;  // lattice before collision (after the last streaming)
    std::vector<int32_t> src_;  // lattice after collision, before streaming
    Well1024a rng_;
    // pow(M, alpha_) for integer M below POW_TABLE_SIZE, filled on first use; -1 = not
    // yet computed (every real value is at least 1). set_options() resets it when alpha_
    // changes.
    // The const collision helpers fill it, so one engine must not be used from two
    // threads at once, not even through its const collision helpers.
    mutable std::vector<double> pow_table_;
};

}  // namespace lgca
