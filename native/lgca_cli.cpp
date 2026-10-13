// native/lgca_cli.cpp — command-line front end of the native engine.
//
// Usage (the reference program's arguments, then optional flags in any position):
//
//   lgca_native <tsteps> <iters> <sens> <dens> [seed] [--W N] [--H N]
//               [--model boson|fermion] [--kernel sum|avg|power] [--alpha A]
//               [--boson-field site|neigh|both] [--boson-align nematic|polar]
//               [--dump-state PATH] [--dump-every N] [--dump-at T1,T2,...]
//               [--precision N]
//
// or, to measure saved lattices with the sweeps' observables (measure_cli.cpp):
//
//   lgca_native measure [--W N --H N] [--factor F] [--snap | --island] [--t T] [--repeat R] FILE...
//
// or, to run or resume one job of the density and alpha sweep as its JavaScript runner
// did, writing the same run.json and state files (run_cli.cpp):
//
//   lgca_native run --size L --alpha A --dens D --seed S --end T DIR
//   lgca_native run --job ID --end T DIR
//
// or, to replay such a job from its seed and measure its island at every sample step
// (trace_cli.cpp):
//
//   lgca_native trace --job ID --end T [--check DIR] [--states DIR] OUT
//
// Without flags it behaves like the reference binary built from lgca/lgca_clean-1.cpp
// in its default configuration (boson model, 120 x 120 lattice): the same arguments,
// parsed the same way (atoi / atof / strtoul; seed from the clock when omitted), and the
// same single output line:
//
//   sens*dens  polar se  nematic se  spatial se  band se  inv_occ se
//
// (each "se" is the standard error over the realizations), averaged at the last time
// step. A sixth positional argument (the reference's wait_ms, used only by its live
// display) is accepted and ignored.
//
// Flags:
//   --W N, --H N        lattice width and height (default 120 x 120; at least 4)
//   --model M           boson (the reference's default build) or fermion
//   --kernel K          sum (default, the reference), avg or power: what the field is
//                       divided by (docs/PHYSICS.md; lgca_engine.hpp, enum Kernel)
//   --alpha A           exponent of the power kernel, a number from 0 to 2 (default 0).
//                       A value other than 0 needs --kernel power: the other kernels
//                       ignore alpha, and a forgotten --kernel would silently run them
//   --boson-field F     site (default, the reference), neigh or both: whose particles
//                       build a boson site's field
//   --boson-align A     nematic (default, the reference) or polar boson kernel
//                       (--boson-field and --boson-align do not apply to the fermion
//                       model, which ignores them, as the Lab does)
//   --dump-state PATH   after the run, write the final state to PATH (format below)
//   --dump-every N      with --dump-state, also write PATH.r<r>.t<t> after every N
//                       steps of every realization r (r from 0, t = steps done)
//   --dump-at LIST      with --dump-state, also write PATH.r<r>.t<t> at each step count
//                       t of the comma-separated LIST; t = 0 is the filled lattice before
//                       the first step
//   --precision N       print the output line with N significant digits instead of the
//                       default 6 (17 digits identify every double exactly, as exact
//                       comparisons with the reference need)
// Flags may also be written --flag=value. Options given more than once: the last wins.
//
// State file: 33 uint32 words of generator state (the 32 state words, then the index,
// as Lgca::get_rng_state() returns them), then the W*H*6 int32 lattice values in the
// engine's layout occ[(j*W + i)*6 + k]. Native byte order, no header.
//
// A file that cannot be written (a dump), or an output line that does not reach stdout
// (stdout closed, read-only, or a pipe whose reader has gone), ends the program with an
// error message and exit status 3, after a requested dump is written.
//
// Inputs on which the reference program has, or can have, undefined behaviour are
// refused here with an error message and exit status 1 (docs/FIDELITY.md): tsteps < 1
// (the reference reads polar[-1]); W or H below 4 (its band scan can read outside the
// lattice; a conservative bound, since an empty lattice never runs the scan); and, when
// a realization runs, a dens the reference cannot convert to an int (not a number, or
// its integer part outside the int range) or, in the boson model, above 357913941, where
// a site's particle count can overflow an int depending on the random fill (also
// conservative; lgca::dens_supported).
//
// The reference measures the observables after every step but prints only those of the
// last step. Measuring draws no random numbers, so this program measures only after the
// last step of each realization. The printed line matches the reference except for
// band, which deliberately counts all six channels (the reference skips channel 0).

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "lgca_engine.hpp"

namespace lgca {
int measure_main(int argc, char** argv);   // measure_cli.cpp
int run_main(int argc, char** argv);       // run_cli.cpp
int trace_main(int argc, char** argv);     // trace_cli.cpp
}

namespace {

template <class T> inline T sqr(T a) { return a * a; }

// Strict integer parsing for the flags (the positional arguments keep the reference's
// lenient atoi / atof on purpose).
bool parse_int(const std::string& text, long min_value, long& value) {
    if (text.empty()) return false;
    char* end = nullptr;
    errno = 0;
    const long v = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || *end != '\0' || v < min_value || v > 1000000000L) return false;
    value = v;
    return true;
}

// Strict parsing of a number for --alpha: the whole text must be a finite number. A
// number too small for a double (strtod reports ERANGE) is accepted as the tiny or zero
// value strtod returns, as JavaScript's Number() reads it; one too large is infinite and
// refused.
bool parse_double(const std::string& text, double& value) {
    if (text.empty()) return false;
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(text.c_str(), &end);
    if ((errno != 0 && errno != ERANGE) || *end != '\0' || !std::isfinite(v)) return false;
    value = v;
    return true;
}

// --dump-at: comma-separated step counts, each from 0 to 10^9, kept sorted.
bool parse_steps(const std::string& text, std::vector<long>& steps) {
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = text.find(',', start);
        long t = 0;
        if (!parse_int(text.substr(start, comma - start), 0, t)) return false;
        steps.push_back(t);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    std::sort(steps.begin(), steps.end());
    return true;
}

// Write generator words and lattice values to a state file (format at the top of this
// file).
bool write_state(const std::string& path, const uint32_t rng[33], const int32_t* occ,
                 std::size_t count) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    bool ok = std::fwrite(rng, sizeof(uint32_t), 33, f) == 33;
    ok = ok && std::fwrite(occ, sizeof(int32_t), count, f) == count;
    ok = (std::fclose(f) == 0) && ok;
    return ok;
}

bool write_state(const std::string& path, const lgca::Lgca& sim) {
    uint32_t rng[33];
    sim.get_rng_state(rng);
    return write_state(path, rng, sim.occ(), sim.cell_count());
}

void usage(const char* program) {
    std::cerr << "usage: " << program
              << " <tsteps> <iters> <sens> <dens> [seed] [wait_ms]"
                 " [--W N] [--H N] [--model boson|fermion]"
                 " [--kernel sum|avg|power] [--alpha A]"
                 " [--boson-field site|neigh|both] [--boson-align nematic|polar]"
                 " [--dump-state PATH] [--dump-every N] [--dump-at T1,T2,...]"
                 " [--precision N]\n";
}

}  // namespace

int main(int argc, char** argv) {
    // A write to a pipe whose reader has gone fails with an error instead of killing the
    // program (SIGPIPE), so every mode reaches its own handling: the simulation still
    // writes a requested dump and exits with status 3, `measure` exits with status 1, and
    // `run` and `trace` finish their files.
    std::signal(SIGPIPE, SIG_IGN);

    // `lgca_native measure ...`: the measuring mode, with its own arguments. (As a
    // simulation's tsteps, "measure" would be refused.)
    if (argc > 1 && std::string(argv[1]) == "measure") return lgca::measure_main(argc - 1, argv + 1);
    if (argc > 1 && std::string(argv[1]) == "run") return lgca::run_main(argc - 1, argv + 1);
    if (argc > 1 && std::string(argv[1]) == "trace") return lgca::trace_main(argc - 1, argv + 1);

    // ---- split the arguments into positionals and --flags --------------------------
    std::vector<std::string> positional;
    lgca::Params params;                 // defaults: 120 x 120, boson
    std::string dump_path;
    long dump_every = 0;
    std::vector<long> dump_at;           // sorted step counts of --dump-at
    long precision = 0;                  // 0 = leave std::cout at its default (6)

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) {
            positional.push_back(arg);
            continue;
        }
        // --name value  or  --name=value
        std::string name = arg.substr(2), value;
        const std::size_t eq = name.find('=');
        if (eq != std::string::npos) {
            value = name.substr(eq + 1);
            name = name.substr(0, eq);
        } else if (a + 1 < argc) {
            value = argv[++a];
        } else {
            std::cerr << "error: --" << name << " needs a value\n";
            return 1;
        }

        long n = 0;
        if (name == "W" || name == "H") {
            if (!parse_int(value, 4, n)) {
                std::cerr << "error: --" << name << " must be an integer of at least 4\n";
                return 1;
            }
            (name == "W" ? params.W : params.H) = static_cast<int>(n);
        } else if (name == "model") {
            if (value == "boson") params.model = lgca::Model::Boson;
            else if (value == "fermion") params.model = lgca::Model::Fermion;
            else {
                std::cerr << "error: --model must be boson or fermion\n";
                return 1;
            }
        } else if (name == "kernel") {
            if (value == "sum") params.kernel = lgca::Kernel::Sum;
            else if (value == "avg") params.kernel = lgca::Kernel::Avg;
            else if (value == "power") params.kernel = lgca::Kernel::Power;
            else {
                std::cerr << "error: --kernel must be sum, avg or power\n";
                return 1;
            }
        } else if (name == "alpha") {
            double alpha = 0;
            if (!parse_double(value, alpha) || alpha < 0 || alpha > 2) {
                std::cerr << "error: --alpha must be a number from 0 to 2\n";
                return 1;
            }
            params.alpha = alpha;
        } else if (name == "boson-field") {
            if (value == "site") params.boson_field = lgca::BosonField::Site;
            else if (value == "neigh") params.boson_field = lgca::BosonField::Neigh;
            else if (value == "both") params.boson_field = lgca::BosonField::Both;
            else {
                std::cerr << "error: --boson-field must be site, neigh or both\n";
                return 1;
            }
        } else if (name == "boson-align") {
            if (value == "nematic") params.boson_align = lgca::BosonAlign::Nematic;
            else if (value == "polar") params.boson_align = lgca::BosonAlign::Polar;
            else {
                std::cerr << "error: --boson-align must be nematic or polar\n";
                return 1;
            }
        } else if (name == "dump-state") {
            dump_path = value;
        } else if (name == "dump-every") {
            if (!parse_int(value, 1, dump_every)) {
                std::cerr << "error: --dump-every must be a positive integer\n";
                return 1;
            }
        } else if (name == "dump-at") {
            dump_at.clear();
            if (!parse_steps(value, dump_at)) {
                std::cerr << "error: --dump-at must be a comma-separated list of step"
                             " counts from 0 to 1000000000\n";
                return 1;
            }
        } else if (name == "precision") {
            if (!parse_int(value, 1, precision) || precision > 40) {
                std::cerr << "error: --precision must be an integer from 1 to 40\n";
                return 1;
            }
        } else {
            std::cerr << "error: unknown option --" << name << "\n";
            usage(argv[0]);
            return 1;
        }
    }
    if (positional.size() < 4) {
        usage(argv[0]);
        return 1;
    }
    if (dump_every > 0 && dump_path.empty()) {
        std::cerr << "error: --dump-every needs --dump-state PATH\n";
        return 1;
    }
    if (!dump_at.empty() && dump_path.empty()) {
        std::cerr << "error: --dump-at needs --dump-state PATH\n";
        return 1;
    }
    if (params.alpha != 0 && params.kernel != lgca::Kernel::Power) {
        std::cerr << "error: --alpha needs --kernel power (the other kernels ignore it)\n";
        return 1;
    }
    if (static_cast<long>(params.W) * params.H > 100000000L) {
        std::cerr << "error: lattice too large (W*H above 10^8)\n";
        return 1;
    }

    // ---- the reference's argument parsing ------------------------------------------
    const int    tsteps = std::atoi(positional[0].c_str());
    const int    iters  = std::atoi(positional[1].c_str());
    const double sens   = std::atof(positional[2].c_str());
    const double dens   = std::atof(positional[3].c_str());
    params.dens = dens;
    params.seed = (positional.size() > 4)
                      ? static_cast<unsigned int>(std::strtoul(positional[4].c_str(), nullptr, 10))
                      : static_cast<unsigned int>(std::time(nullptr));
    if (tsteps < 1) {
        std::cerr << "error: tsteps must be at least 1\n";
        return 1;
    }
    // With iters <= 0 the reference never converts dens, so any dens is fine then.
    if (iters > 0 && !lgca::dens_supported(params.model, dens)) {
        std::cerr << "error: dens must be a number, above -2147483649 and below 2147483648"
                     " (boson model: at most 357913941)\n";
        return 1;
    }

    // ---- realizations ---------------------------------------------------------------
    // The engine is created at the start of the first realization: its constructor
    // seeds the generator and fills the lattice, which the reference does only when a
    // realization starts (with iters <= 0 it fills nothing and draws nothing). Later
    // realizations refill the lattice from the running generator, as the reference does.
    std::optional<lgca::Lgca> sim;

    double polar = 0, polar2 = 0, nemat = 0, nemat2 = 0, entro = 0, entro2 = 0;
    double band = 0, band2 = 0, invocc = 0, invocc2 = 0;

    // Checkpoint files PATH.r<r>.t<t> of --dump-every and --dump-at.
    auto checkpoint = [&](int r) {
        const long t = sim->t();
        const bool wanted = (dump_every > 0 && t > 0 && t % dump_every == 0) ||
                            std::binary_search(dump_at.begin(), dump_at.end(), t);
        if (!wanted) return true;
        const std::string path =
            dump_path + ".r" + std::to_string(r) + ".t" + std::to_string(t);
        if (write_state(path, *sim)) return true;
        std::cerr << "error: cannot write " << path << "\n";
        return false;
    };

    for (int r = 0; r < iters; ++r) {
        if (r == 0) sim.emplace(params);
        else sim->init_lattice();
        if (!checkpoint(r)) return 3;                 // t = 0, only with --dump-at 0
        for (int t = 0; t < tsteps; ++t) {
            sim->step(sens);
            if (!checkpoint(r)) return 3;
        }
        // Observables of the last step, accumulated in the reference's order.
        const lgca::Order o = sim->measure();
        const double b = sim->band();
        polar  += o.polar;    polar2  += sqr(o.polar);
        nemat  += o.nematic;  nemat2  += sqr(o.nematic);
        entro  += o.entropy;  entro2  += sqr(o.entropy);
        band   += b;          band2   += sqr(b);
        invocc += o.inv_occ;  invocc2 += sqr(o.inv_occ);
    }

    // ---- output: the reference's line, the same expressions -----------------------
    auto stderr_of = [iters](double sum, double sum2) {
        return std::sqrt(std::fabs(sqr(sum / iters) - sum2 / iters) / iters);
    };
    const int W = params.W, H = params.H;
    if (precision > 0) std::cout << std::setprecision(static_cast<int>(precision));
    std::cout << sens * dens << " "
              << polar / iters  << " " << stderr_of(polar,  polar2)  << " "
              << nemat / iters  << " " << stderr_of(nemat,  nemat2)  << " "
              << 1.0 - (entro / iters) / std::log(double(W * H)) << " "
              << stderr_of(entro, entro2) / std::log(double(W * H)) << " "
              << band / iters   << " " << stderr_of(band,   band2)   << " "
              << invocc / iters << " " << stderr_of(invocc, invocc2) << "\n";
    std::cout.flush();
    // The line must have reached stdout (a closed or read-only stdout is an error); the
    // final dump is still written.
    const bool printed = static_cast<bool>(std::cout);

    if (!dump_path.empty()) {
        bool ok;
        if (sim) {
            ok = write_state(dump_path, *sim);
        } else {
            // No realization ran (iters <= 0): the reference's state at the end is the
            // freshly seeded generator and its never-filled, all-zero lattice.
            lgca::Well1024a seeded;
            seeded.seed_like_reference(params.seed);
            uint32_t rng[33];
            seeded.get_state(rng);
            const std::vector<int32_t> empty(static_cast<std::size_t>(W) * H * lgca::NODES, 0);
            ok = write_state(dump_path, rng, empty.data(), empty.size());
        }
        if (!ok) {
            std::cerr << "error: cannot write " << dump_path << "\n";
            return 3;
        }
    }
    if (!printed) {
        std::cerr << "error: cannot write the output line\n";
        return 3;
    }
    return 0;
}
