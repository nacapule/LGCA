// native/wasm/lgca_wasm.cpp — a plain C interface to the engine (native/lgca_engine.hpp)
// for the WebAssembly build (native/build-wasm.sh).
//
// The engine is the same C++ the native program runs, compiled unchanged; this file only
// adds functions a JavaScript page can call. One module instance holds at most one engine
// (a second lgca_create() replaces the first). Several engines need several module
// instances, which cost nothing extra to load from the same glue.
//
// Every number crosses the boundary as a plain int, unsigned int or double, and arrays
// cross as pointers into the module's memory (the caller allocates them with _malloc, or
// reads the engine's own lattice in place through lgca_occ_ptr()). From JavaScript:
//
//   const M = await createLgcaWasm();                      // the glue's factory
//   M._lgca_create(W, H, model, dens, seed, kernel, alpha, field, align) === 1
//   M._lgca_step(n, sens);                                 // n steps at sensitivity sens
//   M._lgca_set_options(kernel, alpha, field, align);      // switch options between steps
//   new Int32Array(M.HEAP32.buffer, M._lgca_occ_ptr(), M._lgca_cell_count())
//
// Option codes (the order of the enums in lgca_engine.hpp; the first value of each is
// the reference program):
//   model   0 fermion, 1 boson
//   kernel  0 sum, 1 avg, 2 power (exponent alpha, 0..2)
//   field   0 site, 1 neigh, 2 both        (boson only)
//   align   0 nematic, 1 polar             (boson only)
//
// No function here throws or aborts on bad input. lgca_create() checks the parameters
// first (with the engine's own check, lgca::params_error) and returns 0 with a message
// for lgca_error(); the other functions do nothing (or return 0) when no engine exists.
// The build has no C++ exception catching, so a throw would end the module: the checks
// make sure the engine is never asked to throw. lgca_create() also checks that the
// memory for the engine can be had before it builds it (see there).
//
// Layouts, identical to the native program's state files and to the Lab's engine:
//   lattice  int32, occ[(j*W + i)*6 + k] = particles at site (i, j) in channel k
//   rng      33 uint32 words: the 32 WELL1024a state words, then the index (0..31)
//   measure  5 doubles: polar, nematic, entropy, 1/occupied (0 if none), occupied

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>

#include <emscripten/emscripten.h>

#include "../lgca_engine.hpp"

namespace {

std::unique_ptr<lgca::Lgca> engine;
const char* last_error = "";

// The lattice must fit comfortably in WebAssembly memory (at most 2 GiB with the build's
// settings; the engine holds two int32 lattices, occ and src). 2^27 cells = 512 MiB per
// lattice, for example 4700 x 4700 sites. This also keeps W*H*6 far from overflowing the
// 32-bit size_t of WebAssembly.
constexpr double MAX_CELLS = 134217728.0;   // 2^27

// Whether malloc can supply the engine's memory right now: the engine object, its table of
// powers (lgca::POW_TABLE_SIZE doubles) and its two lattices of `cells` int32 values, in
// the order the engine allocates them. Allocates the same blocks, then frees them. With
// memory growth on, malloc returns null when the memory cannot grow (the browser refuses
// it, or the 2 GiB ceiling is reached), whereas the engine's own allocations (operator new,
// std::vector) would abort the module, as the build has no exception catching. A check,
// not a guarantee: it cannot see another allocation that happens in between, but nothing
// else runs in this module between the check and the engine's construction. The pointers
// are volatile because the compiler may otherwise remove a malloc whose block is never
// used, and with it the check. After construction the engine allocates nothing more:
// stepping and option changes (lgca_set_options) reuse these blocks.
bool memory_available(std::size_t cells) {
    void* volatile object = std::malloc(sizeof(lgca::Lgca));
    void* volatile powers = std::malloc(lgca::POW_TABLE_SIZE * sizeof(double));
    void* volatile occ = std::malloc(cells * sizeof(int32_t));
    void* volatile src = std::malloc(cells * sizeof(int32_t));
    const bool ok = object && powers && occ && src;
    std::free(src);
    std::free(occ);
    std::free(powers);
    std::free(object);
    return ok;
}

// N and t are int64 in the engine and cross the boundary as doubles, exact for integers
// up to 2^53 in size.
constexpr double MAX_EXACT = 9007199254740992.0;   // 2^53
bool exact_integer(double x) { return std::isfinite(x) && x == std::trunc(x) && std::fabs(x) <= MAX_EXACT; }

}  // namespace

extern "C" {

// Build an engine: seed the generator and fill the lattice for the first realization,
// exactly as lgca::Lgca's constructor does. Returns 1, or 0 if the parameters are refused
// or the memory for the lattices cannot be had (any earlier engine is then gone too;
// lgca_error() says why).
EMSCRIPTEN_KEEPALIVE
int lgca_create(int W, int H, int model, double dens, uint32_t seed,
                int kernel, double alpha, int field, int align) {
    engine.reset();
    last_error = "";
    if (model < 0 || model > 1) { last_error = "model must be 0 (fermion) or 1 (boson)"; return 0; }
    if (kernel < 0 || kernel > 2) { last_error = "kernel must be 0 (sum), 1 (avg) or 2 (power)"; return 0; }
    if (field < 0 || field > 2) { last_error = "field must be 0 (site), 1 (neigh) or 2 (both)"; return 0; }
    if (align < 0 || align > 1) { last_error = "align must be 0 (nematic) or 1 (polar)"; return 0; }

    lgca::Params p;
    p.W = W;
    p.H = H;
    p.model = model == 0 ? lgca::Model::Fermion : lgca::Model::Boson;
    p.dens = dens;
    p.seed = seed;
    p.kernel = kernel == 0 ? lgca::Kernel::Sum : kernel == 1 ? lgca::Kernel::Avg : lgca::Kernel::Power;
    p.alpha = alpha;
    p.boson_field = field == 0 ? lgca::BosonField::Site
                  : field == 1 ? lgca::BosonField::Neigh : lgca::BosonField::Both;
    p.boson_align = align == 0 ? lgca::BosonAlign::Nematic : lgca::BosonAlign::Polar;

    if (const char* error = lgca::params_error(p)) { last_error = error; return 0; }
    if (double(W) * double(H) * lgca::NODES > MAX_CELLS) {
        last_error = "lattice too large for WebAssembly memory (W*H*6 must be at most 2^27)";
        return 0;
    }
    if (!memory_available(static_cast<std::size_t>(W) * H * lgca::NODES)) {
        last_error = "not enough memory for this lattice (the WebAssembly memory could not grow)";
        return 0;
    }
    engine = std::make_unique<lgca::Lgca>(p);
    return 1;
}

// Free the engine and its lattices.
EMSCRIPTEN_KEEPALIVE
void lgca_destroy(void) { engine.reset(); }

// Why the last lgca_create() refused its parameters ("" if it did not). A C string in
// the module's memory; read it with UTF8ToString().
EMSCRIPTEN_KEEPALIVE
const char* lgca_error(void) { return last_error; }

// Start a new realization: refill the lattice from the current generator state (the
// reference does this at the start of every realization, without reseeding).
EMSCRIPTEN_KEEPALIVE
void lgca_init_lattice(void) {
    if (engine) engine->init_lattice();
}

// n time steps at sensitivity sens (each: collision at every occupied site, then
// streaming). Looping here rather than in JavaScript saves one call per step.
EMSCRIPTEN_KEEPALIVE
void lgca_step(int n, double sens) {
    if (!engine) return;
    for (int s = 0; s < n; ++s) engine->step(sens);
}

// The lattice after the last step (occ) and the collision output that step streamed
// (src), W*H*6 int32 values each. The pointers stay valid until lgca_create() or
// lgca_destroy(); a JavaScript view on them must be rebuilt after the memory grows.
EMSCRIPTEN_KEEPALIVE
int32_t* lgca_occ_ptr(void) { return engine ? engine->occ() : nullptr; }

EMSCRIPTEN_KEEPALIVE
const int32_t* lgca_src_ptr(void) { return engine ? engine->src() : nullptr; }

// Number of int32 values in each lattice: W*H*6.
EMSCRIPTEN_KEEPALIVE
int lgca_cell_count(void) { return engine ? static_cast<int>(engine->cell_count()) : 0; }

// Replace the lattice with the W*H*6 values at `in`. To replay a saved state, also set
// the generator (lgca_set_rng), N (lgca_set_nuparts) and the step counter (lgca_set_t) to
// the saved ones: the engine then continues exactly as from that state, and measures as
// it would have. N matters even when the state comes from the same run with the same
// settings: it depends on the seed and the realization (the random fill adds particles),
// and it divides the observables. The collision output (src) is not part of a state: the
// next step rewrites it entirely.
EMSCRIPTEN_KEEPALIVE
void lgca_set_occ(const int32_t* in) {
    if (engine && in) std::memcpy(engine->occ(), in, engine->cell_count() * sizeof(int32_t));
}

// The generator state, 33 words (32 state words, then the index), out of or into `ptr`.
EMSCRIPTEN_KEEPALIVE
void lgca_get_rng(uint32_t* out) {
    if (engine && out) engine->get_rng_state(out);
}

EMSCRIPTEN_KEEPALIVE
void lgca_set_rng(const uint32_t* in) {
    if (engine && in) engine->set_rng_state(in);
}

// The reference's measure_order() on the current lattice: writes 5 doubles to `out`
// (polar, nematic, entropy, 1/occupied, occupied). Draws no random numbers.
EMSCRIPTEN_KEEPALIVE
void lgca_measure(double* out) {
    if (!engine || !out) return;
    const lgca::Order m = engine->measure();
    out[0] = m.polar;
    out[1] = m.nematic;
    out[2] = m.entropy;
    out[3] = m.inv_occ;
    out[4] = static_cast<double>(m.occupied);
}

// The reference's measure_band() on the current lattice. Draws no random numbers.
EMSCRIPTEN_KEEPALIVE
double lgca_band(void) { return engine ? engine->band() : 0.0; }

// N, the particle count the observables divide by, and the steps since the last
// lattice fill. Doubles, exact up to 2^53.
EMSCRIPTEN_KEEPALIVE
double lgca_nuparts(void) { return engine ? static_cast<double>(engine->nuparts()) : 0.0; }

EMSCRIPTEN_KEEPALIVE
double lgca_t(void) { return engine ? static_cast<double>(engine->t()) : 0.0; }

// Set N or the step counter, for replaying a saved state (see lgca_set_occ). Neither
// enters the dynamics. Returns 1, or 0 (and changes nothing) without an engine or if the
// value is not an integer of size at most 2^53, or for t, if it is negative. N may be
// negative: the reference counts a negative dens that way (int(dens) * 6 * W * H).
EMSCRIPTEN_KEEPALIVE
int lgca_set_nuparts(double n) {
    if (!engine || !exact_integer(n)) return 0;
    engine->set_nuparts(static_cast<std::int64_t>(n));
    return 1;
}

EMSCRIPTEN_KEEPALIVE
int lgca_set_t(double t) {
    if (!engine || !exact_integer(t) || t < 0) return 0;
    engine->set_t(static_cast<std::int64_t>(t));
    return 1;
}

// Change the research options of the current engine between steps (codes as for
// lgca_create; the Lab does this when its user switches an option during a run, and the
// change acts from the next step on). Returns 1, or 0 (and changes nothing) without an
// engine, for an unknown code, or if alpha is not a number from 0 to 2.
EMSCRIPTEN_KEEPALIVE
int lgca_set_options(int kernel, double alpha, int field, int align) {
    if (!engine) return 0;
    if (kernel < 0 || kernel > 2 || field < 0 || field > 2 || align < 0 || align > 1) return 0;
    if (!(alpha >= 0.0 && alpha <= 2.0)) return 0;   // set_options() would throw; also NaN
    engine->set_options(kernel == 0 ? lgca::Kernel::Sum : kernel == 1 ? lgca::Kernel::Avg : lgca::Kernel::Power,
                        alpha,
                        field == 0 ? lgca::BosonField::Site
                        : field == 1 ? lgca::BosonField::Neigh : lgca::BosonField::Both,
                        align == 0 ? lgca::BosonAlign::Nematic : lgca::BosonAlign::Polar);
    return 1;
}

}  // extern "C"
