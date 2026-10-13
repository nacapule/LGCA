// native/run_cli.cpp — the native program's `run` mode: one sweep job, run as the
// density and alpha sweep's JavaScript runner ran it.
//
// Usage:
//
//   lgca_native run --size L --alpha A --dens D --seed S --end T DIR
//   lgca_native run --job ID --end T DIR
//
// It does what that runner did for one job (`job ID T`), and writes the same files:
//
//   DIR/run.json            the job, the `series` rows (one per sample time) and the
//                           `snaps` (one per checkpoint, thresholdChecks included);
//   DIR/state-t<T>.bin.gz   the lattice at each checkpoint T: the L*L*6 int32 values,
//                           little-endian, gzip-compressed. Its `snaps` entry holds the
//                           SHA-256 of the uncompressed bytes and the 33 generator words.
//
// The job's fixed settings are those of the sweep: boson model, power kernel with
// exponent alpha, field from the site and its 6 neighbours ("both"), polar alignment,
// sensitivity 6, an L x L lattice. ID is the sweep's job name,
// both-L<L>-a<A>-d<D>-seed<S> (both-L90-a0.85-d0.4-seed12345); with --size etc. the
// program builds that name itself. T must be one of the checkpoints:
//
//   0 100 500 1000 2000 5000 10000 20000 40000 60000 80000, then every 20000 up to 2560000
//
// Sample times (a `series` row): every 25 steps up to step 5000, then every 100.
//
// Resume: if DIR/run.json exists, the job continues from its last checkpoint, as the
// JavaScript runner does. The program checks that the file belongs to the same job and
// engine, that its `snaps` are every checkpoint up to the last one and its `series` rows
// up to that checkpoint every sample step, each once and in order, every one of them with
// its measurements (the 18 series fields, finite numbers, whole ones where they count
// something; a snap also its thresholdChecks, the lattice's SHA-256 and the 33 generator
// words), so that a folder counted as done holds all its measurements. It loads that
// checkpoint's lattice (refusing it if its SHA-256 differs from the saved one) and
// generator words, drops the `series` rows written after that checkpoint (they are
// measured again), and continues to T. If the last checkpoint is already at or past T
// there is nothing to do: the program says so and leaves the folder as it is.
//
// The written run.json is the JavaScript runner's, with one more field, "engine", right
// after "engineHash": which native build wrote the file (the SHA-256 of the native/ sources,
// the compiler and its flags). "engineHash" keeps the pinned value of the JavaScript
// engine that the runner checks, because this program reproduces that engine's dynamics
// exactly; a folder written by either program can therefore be continued by the other (by
// the runner only up to step 80000, its last checkpoint).
// The numbers are written as JSON.stringify writes them; rows and snaps kept from an
// existing run.json are copied as they are written there.
//
// A folder without run.json is a new job, unless it holds state files: then run.json,
// which says whose they are, is missing, and the program refuses the folder rather than
// write over them. The one exception is a lone state-t0.bin.gz holding this job's
// starting lattice, byte for byte: a first save stopped before its run.json (the job
// starts again, writing the same file).
//
// One program at a time: the program locks DIR itself (flock on the folder, so the
// folder gets no extra file) before it reads anything there, and holds the lock until it
// ends. A second `lgca_native run` on the same folder, or an `lgca_native trace` that
// writes there, stops at once with an error. (The JavaScript runner takes no lock: do
// not run it on a folder a native program is working in.)
//
// Files are replaced atomically (written to a temporary file of their own, NAME.tmp.
// <process id>.<n>, flushed to the disk, renamed), and run.json is written after the
// state file it points to, so the program can be stopped at any moment: the folder then
// always holds a complete run.json whose last checkpoint is on disk. A stop during a
// write can leave a stray temporary file; the next run that goes on to write in the
// folder removes those of run.json and the state files first (also NAME.tmp, the name
// earlier builds and the JavaScript runner use).
//
// Prints "<ID>: t=<T>, <seconds>s" at the end, as the JavaScript runner does. Errors stop
// the program with a message and exit status 1. Options may also be written --option=value.

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "file_io.hpp"
#include "json.hpp"
#include "lgca_engine.hpp"
#include "measure.hpp"
#include "sha256.hpp"
#include "sweep_job.hpp"

namespace lgca {

namespace {

// The run.json text, in the JavaScript runner's key order:
// {job:{...job,end}, engineHash, engine, model, alignment, stage, seconds, series, snaps}.
// `series` and `snaps` are the array contents, comma-separated, without the brackets.
std::string run_json(const Job& job, long end, double seconds, const std::string& series,
                     const std::string& snaps) {
    std::string s;
    s.reserve(series.size() + snaps.size() + 1024);
    s += "{\"job\":{\"id\":" + json_string(job.id);
    s += ",\"size\":" + json_number(static_cast<double>(job.size));
    s += ",\"alpha\":" + json_number(job.alpha);
    s += ",\"dens\":" + json_number(job.dens);
    s += ",\"seed\":" + json_number(job.seed);
    s += ",\"sens\":" + json_number(SENS);
    s += ",\"field\":\"both\"";
    s += ",\"end\":" + json_number(static_cast<double>(end));
    s += "},\"engineHash\":" + json_string(ENGINE_HASH);
    s += ",\"engine\":" + engine_json("lgca_native run");
    s += ",\"model\":\"boson\",\"alignment\":\"polar\",\"stage\":\"after streaming\"";
    s += ",\"seconds\":" + json_number(seconds);
    s += ",\"series\":[" + series + "],\"snaps\":[" + snaps + "]}";
    return s;
}


// The measurements a series row and a snap hold (series_json): each a finite number, a
// whole one from 0 where it counts something.
struct Field {
    const char* name;
    bool count;
};
const Field SERIES_FIELDS[] = {
    {"N", true},           {"mean", false},          {"nmax", true},
    {"kmax", true},        {"siteFraction", false},  {"channelFraction", false},
    {"effectiveArea", false}, {"occupied", true},    {"polar", false},
    {"localPolar", false}, {"weightedLocalPolar", false}, {"threshold", true},
    {"clusters", true},    {"clusterMass", true},    {"clusterArea", true},
    {"clusterDensity", false}, {"clusterFraction", false}, {"clusterPolar", false}};
static_assert(sizeof SERIES_FIELDS / sizeof SERIES_FIELDS[0] == 18, "the 18 fields of series_json");
// Counts can pass 2^53 (a sum over 10^8 sites of up to 2^31 - 1 particles): every finite
// whole number is allowed.
constexpr double MAX_COUNT = std::numeric_limits<double>::max();

bool is_finite_number(const JsonValue* v) { return v && v->is_number() && std::isfinite(v->number); }

// The first measurement `entry` lacks or holds as something else, or "" if they are all
// there; with `snap`, also its thresholdChecks (factors 1.5 and 4, snap_json).
std::string bad_measurement(const JsonValue& entry, bool snap) {
    for (const Field& f : SERIES_FIELDS) {
        const JsonValue* v = entry.get(f.name);
        if (f.count ? !is_integer(v, 0, MAX_COUNT) : !is_finite_number(v)) return f.name;
    }
    if (!snap) return "";
    const JsonValue* checks = entry.get("thresholdChecks");
    if (!checks || !checks->is_array() || checks->items.size() != 2) return "thresholdChecks";
    for (std::size_t c = 0; c < 2; ++c) {
        const JsonValue& check = checks->items[c];
        const JsonValue* factor = check.get("factor");
        if (!is_finite_number(factor) || factor->number != THRESHOLD_CHECK_FACTORS[c] ||
            !is_integer(check.get("threshold"), 0, MAX_COUNT) || !is_integer(check.get("mass"), 0, MAX_COUNT) ||
            !is_integer(check.get("area"), 0, MAX_COUNT) || !is_finite_number(check.get("density")))
            return "thresholdChecks";
    }
    return "";
}

const char* const MEASUREMENTS_TEXT =
    " (each measurement a finite number, a whole one where it counts something)";

// Where a resumed job continues.
struct Checkpoint {
    long t = -1;                 // the last checkpoint's step; -1: no run.json, start at 0
    double previous_seconds = 0; // run.json's "seconds"
    std::string series;          // the kept series rows (t <= the checkpoint), as written
    std::string snaps;           // every snaps entry, as written
};

// Read DIR/run.json and restore the job's last checkpoint into `sim` (as the JavaScript
// runner resumes a job, with stricter checks of the file's contents). Throws
// std::runtime_error with the reason if the file cannot be used.
Checkpoint resume(const std::string& dir, const Job& job, Lgca& sim) {
    namespace fs = std::filesystem;
    Checkpoint cp;
    const std::string path = (fs::path(dir) / "run.json").string();
    std::vector<unsigned char> bytes;
    std::string error;
    // A run.json of a job to step 2560000, the last checkpoint, holds 25751 series rows and
    // 135 snaps, about 12 MB at any lattice size; the limit keeps a damaged file from
    // filling memory.
    if (!read_file(path, std::size_t(64) << 20, bytes, error))
        throw std::runtime_error(path + ": cannot read: " + error);
    const std::string text(bytes.begin(), bytes.end());
    JsonValue old;
    if (!parse_json(text, old, error)) throw std::runtime_error(path + ": not valid JSON: " + error);
    const auto fail = [&path](const std::string& why) { return std::runtime_error(path + ": " + why); };
    if (!old.is_object()) throw fail("not a JSON object");

    // The same engine and the same job (the runner's "Engine changed", "Settings changed").
    const JsonValue* engine_hash = old.get("engineHash");
    if (!engine_hash || !engine_hash->is_string() || engine_hash->text != ENGINE_HASH)
        throw fail("Engine changed: engineHash is not the pinned engine's");
    const JsonValue* old_job = old.get("job");
    if (!old_job || !old_job->is_object()) throw fail("no \"job\" object");
    const JsonValue* old_id = old_job->get("id");
    if (!old_id || !old_id->is_string() || old_id->text != job.id)
        throw fail("Settings changed: the saved job's \"id\" is not " + job.id);
    const auto same_number = [old_job](const char* key, double value) {
        const JsonValue* v = old_job->get(key);
        return v && v->is_number() && v->number == value;
    };
    const JsonValue* field = old_job->get("field");
    if (!same_number("size", static_cast<double>(job.size)) || !same_number("alpha", job.alpha) ||
        !same_number("dens", job.dens) || !same_number("seed", job.seed) ||
        !same_number("sens", SENS) || !field || !field->is_string() || field->text != "both")
        throw fail("Settings changed: the saved job is not " + job.id);
    const JsonValue* model = old.get("model");
    const JsonValue* alignment = old.get("alignment");
    if (!model || !model->is_string() || model->text != "boson" || !alignment ||
        !alignment->is_string() || alignment->text != "polar")
        throw fail("Settings changed: model and alignment must be boson and polar");

    // The snaps: checkpoints in increasing order; the last one is where the job continues.
    const JsonValue* snaps = old.get("snaps");
    if (!snaps || !snaps->is_array() || snaps->items.empty())
        throw fail("no checkpoint (\"snaps\" is missing or empty)");
    // Every snap is a checkpoint step (the runner writes snaps only there) with its file,
    // the lattice's SHA-256, the 33 generator words and its measurements.
    long previous_t = -1;
    uint32_t words[33];
    for (const JsonValue& snap : snaps->items) {
        const JsonValue* t = snap.get("t");
        if (!snap.is_object() || !is_integer(t, 0, 1e9) || !is_tick(static_cast<long>(t->number)) ||
            !(t->number > previous_t))
            throw fail(std::string("\"snaps\" entries must be objects whose \"t\" are increasing"
                                   " checkpoint steps (") + CHECKPOINT_TEXT + ")");
        previous_t = static_cast<long>(t->number);
        const std::string at = "the checkpoint at t=" + std::to_string(previous_t);
        const std::string name = "state-t" + std::to_string(previous_t) + ".bin.gz";
        const JsonValue* file = snap.get("file");
        if (!file || !file->is_string() || file->text != name)
            throw fail(at + " must have \"file\" " + name);
        const JsonValue* sha = snap.get("sha256");
        if (!sha || !sha->is_string() || !is_lower_hex(sha->text, 64))
            throw fail(at + ": its \"sha256\" is not 64 lowercase hex digits");
        const JsonValue* rng = snap.get("rng");
        if (!rng || !rng->is_array() || rng->items.size() != 33)
            throw fail(at + ": its \"rng\" is not 33 numbers");
        for (int i = 0; i < 33; ++i) {
            // 32 state words, then the index, which is below 32.
            if (!is_integer(&rng->items[i], 0, i < 32 ? 4294967295.0 : 31.0))
                throw fail(at + ": its \"rng\" holds an invalid generator word");
            words[i] = static_cast<uint32_t>(rng->items[i].number);
        }
        const std::string bad = bad_measurement(snap, true);
        if (!bad.empty()) throw fail(at + " has no valid \"" + bad + "\"" + MEASUREMENTS_TEXT);
        append_item(cp.snaps, text.substr(snap.begin, snap.end - snap.begin));
    }
    // Every checkpoint up to the last one, none skipped (the loop above checked the order).
    const std::vector<long> ticks = tick_steps(previous_t);
    if (snaps->items.size() != ticks.size())
        throw fail("\"snaps\" must hold every checkpoint step from 0 to the last one, t=" +
                   std::to_string(previous_t) + " (" + std::to_string(ticks.size()) + " of them), not " +
                   std::to_string(snaps->items.size()));
    // The last checkpoint, where the job continues (`words` holds its generator words).
    const JsonValue& last = snaps->items.back();
    cp.t = previous_t;
    const std::string file_name = "state-t" + std::to_string(cp.t) + ".bin.gz";
    const JsonValue* sha = last.get("sha256");
    const JsonValue* last_n = last.get("N");

    // The lattice: its SHA-256 must be the saved one (the runner's "State hash mismatch").
    const std::string state_path = (fs::path(dir) / file_name).string();
    const std::size_t expected = sim.cell_count() * sizeof(int32_t);
    if (!read_file(state_path, expected, bytes, error))
        throw std::runtime_error(state_path + ": cannot read: " + error);
    if (bytes.size() != expected)
        throw std::runtime_error(state_path + ": " + std::to_string(bytes.size()) +
                                 " bytes; an L = " + std::to_string(job.size) + " lattice has " +
                                 std::to_string(expected));
    if (sha256_hex(bytes.data(), bytes.size()) != sha->text)
        throw std::runtime_error(state_path + ": State hash mismatch: its SHA-256 is not the one"
                                              " run.json saved");
    std::memcpy(sim.occ(), bytes.data(), expected);
    sim.set_rng_state(words);
    // The particle count of a freshly filled lattice with this seed; particles are
    // conserved, so the checkpoint must have the same (the runner's "Particle count
    // mismatch").
    if (last_n->number != static_cast<double>(sim.nuparts()))
        throw fail("Particle count mismatch: the last checkpoint has N = " +
                   json_number(last_n->number) + ", the job has " + std::to_string(sim.nuparts()));
    // The loaded lattice itself must hold that many particles, none negative, and no site
    // more than a 32-bit integer holds (the engine adds a site's channels in one).
    const int32_t* occ = sim.occ();
    int64_t total = 0;
    for (std::size_t x = 0; x < sim.cell_count() / 6; ++x) {
        int64_t site = 0;
        for (int k = 0; k < 6; ++k) {
            if (occ[6 * x + k] < 0) throw std::runtime_error(state_path + ": a negative count");
            site += occ[6 * x + k];
        }
        if (site > INT32_MAX)
            throw std::runtime_error(state_path + ": a site holds more than 2^31 - 1 particles");
        total += site;
    }
    if (total != static_cast<int64_t>(sim.nuparts()))
        throw std::runtime_error(state_path + ": Particle count mismatch: the lattice holds " +
                                 std::to_string(total) + ", the job has " +
                                 std::to_string(sim.nuparts()));

    // The series rows up to the checkpoint; later rows are measured again (the runner's
    // series = old.series.filter(x => x.t <= last.t)).
    // Those kept rows must be every sample step up to the checkpoint, each once, in order,
    // each with its measurements.
    const JsonValue* series = old.get("series");
    if (!series || !series->is_array()) throw fail("no \"series\" array");
    const std::vector<long> samples = sample_steps(cp.t);
    std::size_t kept = 0;
    for (const JsonValue& row : series->items) {
        const JsonValue* t = row.get("t");
        if (!row.is_object() || !t || !t->is_number())
            throw fail("\"series\" rows must be objects with a number \"t\"");
        if (t->number > cp.t) continue;
        if (kept >= samples.size() || t->number != static_cast<double>(samples[kept]))
            throw fail("the \"series\" rows up to the last checkpoint, t=" + std::to_string(cp.t) +
                       ", must be every sample step (every 25 steps to 5000, then every 100), each"
                       " once and in order: found t=" + json_number(t->number) + " where " +
                       (kept < samples.size() ? "t=" + std::to_string(samples[kept]) : std::string("none")) +
                       " belongs");
        const std::string bad = bad_measurement(row, false);
        if (!bad.empty())
            throw fail("the \"series\" row at t=" + std::to_string(samples[kept]) + " has no valid \"" + bad +
                       "\"" + MEASUREMENTS_TEXT);
        ++kept;
        append_item(cp.series, text.substr(row.begin, row.end - row.begin));
    }
    if (kept != samples.size())
        throw fail("the \"series\" rows up to the last checkpoint, t=" + std::to_string(cp.t) +
                   ", must be every sample step (every 25 steps to 5000, then every 100): the rows"
                   " stop before t=" + std::to_string(samples[kept]));

    // The runner's previousSeconds = old.seconds || 0.
    const JsonValue* seconds = old.get("seconds");
    if (seconds && seconds->is_number() && std::isfinite(seconds->number))
        cp.previous_seconds = seconds->number;
    return cp;
}

// A folder without run.json (see the top of this file): refused if it holds state files,
// unless the only one is state-t0.bin.gz holding this job's starting lattice. `sim` is
// the job at t = 0. Throws std::runtime_error with the reason.
void check_new_folder(const std::string& dir, const Lgca& sim) {
    namespace fs = std::filesystem;
    std::vector<std::string> states;
    for (const std::string& name : folder_entries(dir))
        if (is_state_file_name(name)) states.push_back(name);
    if (states.empty()) return;
    if (states.size() == 1 && states[0] == "state-t0.bin.gz") {
        const std::string path = (fs::path(dir) / states[0]).string();
        const std::size_t expected = sim.cell_count() * sizeof(int32_t);
        std::vector<unsigned char> bytes;
        std::string error;
        std::string why;
        if (!read_file(path, expected, bytes, error))
            why = "cannot be read (" + error + ")";
        else if (bytes.size() != expected)
            why = "holds " + std::to_string(bytes.size()) + " bytes, not the " +
                  std::to_string(expected) + " of this job's lattice";
        else if (std::memcmp(bytes.data(), sim.occ(), expected) != 0)
            why = "is not the starting lattice of this job";
        else
            return;   // a first save stopped before its run.json: the job starts again
        throw std::runtime_error(path + ": there is no run.json, and this file " + why +
                                 ": refusing to write over it. Put the folder's run.json back,"
                                 " or move the file away to start the job afresh");
    }
    std::string list;
    for (const std::string& name : states) list += (list.empty() ? "" : " ") + name;
    throw std::runtime_error(dir + ": holds saved lattices (" + list + ") but no run.json, which"
                             " says whose they are: refusing to start the job over them. Put"
                             " run.json back to continue the job, or move the state files away"
                             " to start it afresh");
}

// Remove the stray temporary files of run.json and of the state files that a stopped
// run left in `dir` (atomic_temp_target). The caller holds the folder's lock, so no other
// run is writing them. Throws std::runtime_error if one cannot be removed.
void remove_stray_temp_files(const std::string& dir) {
    namespace fs = std::filesystem;
    for (const std::string& name : folder_entries(dir)) {
        const std::string target = atomic_temp_target(name);
        if (target != "run.json" && !is_state_file_name(target)) continue;
        const fs::path path = fs::path(dir) / name;
        std::error_code ec;
        if (fs::is_directory(fs::symlink_status(path, ec))) continue;
        if (!fs::remove(path, ec) && ec)
            throw std::runtime_error(path.string() + ": cannot remove this stray temporary file: " +
                                     ec.message());
    }
}

void usage() {
    std::cerr << "usage: lgca_native run --size L --alpha A --dens D --seed S --end T DIR\n"
                 "       lgca_native run --job ID --end T DIR\n";
}

// Run the job (the body of run_main, which turns exceptions into error messages).
int run_job(const Job& job, long end, const std::string& dir) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec || !fs::is_directory(dir)) {
        std::cerr << "error: " << dir << ": cannot create the folder"
                  << (ec ? ": " + ec.message() : std::string()) << "\n";
        return 1;
    }
    // Before anything in the folder is read.
    const FolderLock lock(dir, dir + ": another program (lgca_native run or trace) is working in"
                                     " this folder; run this again when it has finished");

    const Params params = job_params(job);
    Lgca sim(params);   // seeded and filled, as new LGCA(...) in the runner
    Measurer measurer(params.W, params.H);
    const std::string run_path = (fs::path(dir) / "run.json").string();

    // run.json counts as there when its name is, even as a link to nowhere (then it
    // cannot be read, and the folder is refused rather than started afresh).
    Checkpoint cp;
    const fs::file_status run_status = fs::symlink_status(run_path, ec);
    if (ec && run_status.type() != fs::file_type::not_found)
        throw std::runtime_error(run_path + ": " + ec.message());
    if (fs::exists(run_status)) cp = resume(dir, job, sim);
    else check_new_folder(dir, sim);
    if (cp.t >= end) {
        std::cout << job.id << ": already at t=" << cp.t << ", nothing to do\n";
        return 0;
    }
    // About to write: first clear what a stopped run left.
    remove_stray_temp_files(dir);

    const auto begin = std::chrono::steady_clock::now();
    // The runner's previousSeconds + (Date.now() - begin) / 1000, with whole milliseconds.
    const auto seconds = [&begin, &cp]() {
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - begin).count();
        return cp.previous_seconds + static_cast<double>(ms) / 1000;
    };
    const std::size_t cells = sim.cell_count();
    const std::size_t lattice_bytes = cells * sizeof(int32_t);
    std::vector<unsigned char> gz;
    std::string error;

    for (long t = cp.t + 1; t <= end; ++t) {
        if (t) sim.step(SENS);
        if (!is_sample(t)) continue;
        const int32_t* occ = sim.occ();
        const Measurements m = measurer.measure(occ, 2);
        bool negative = false;
        for (std::size_t c = 0; c < cells; ++c) negative |= occ[c] < 0;
        if (m.N != static_cast<double>(sim.nuparts()) || negative) {
            std::cerr << "error: " << job.id << ": Invalid particle count at t=" << t << "\n";
            return 1;
        }
        if (is_tick(t)) {
            // The lattice's raw bytes: its SHA-256, and the gzip file.
            const std::string file = "state-t" + std::to_string(t) + ".bin.gz";
            const std::string sha = sha256_hex(occ, lattice_bytes);
            if (!gzip_compress(occ, lattice_bytes, gz, error) ||
                !write_file_atomic((fs::path(dir) / file).string(), gz.data(), gz.size(), error)) {
                std::cerr << "error: " << error << "\n";
                return 1;
            }
            uint32_t rng[33];
            sim.get_rng_state(rng);
            append_item(cp.snaps, snap_json(t, file, sha, rng, m, measurer.threshold_checks(occ)));
        }
        append_item(cp.series, series_json(t, m));
        if (is_tick(t)) {
            const std::string text = run_json(job, end, seconds(), cp.series, cp.snaps);
            if (!write_file_atomic(run_path, text.data(), text.size(), error)) {
                std::cerr << "error: " << error << "\n";
                return 1;
            }
        }
    }
    char line[64];
    std::snprintf(line, sizeof line, "%.1f", seconds());
    std::cout << job.id << ": t=" << end << ", " << line << "s\n";
    return 0;
}

}  // namespace

int run_main(int argc, char** argv) {
    std::string job_name, dir;
    bool have_size = false, have_alpha = false, have_dens = false, have_seed = false;
    long long size = 0, seed = 0, end = -1;
    double alpha = 0, dens = 0;

    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) {
            if (!dir.empty()) {
                std::cerr << "error: give one folder (also got " << arg << ")\n";
                usage();
                return 1;
            }
            dir = arg;
            continue;
        }
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
        if (name == "job") {
            job_name = value;
        } else if (name == "size") {
            if (!parse_integer(value, 4, 10000, size)) {
                std::cerr << "error: --size must be an integer from 4 to 10000\n";
                return 1;
            }
            have_size = true;
        } else if (name == "alpha") {
            if (!parse_number(value, alpha)) {
                std::cerr << "error: --alpha must be a number\n";
                return 1;
            }
            have_alpha = true;
        } else if (name == "dens") {
            if (!parse_number(value, dens)) {
                std::cerr << "error: --dens must be a number\n";
                return 1;
            }
            have_dens = true;
        } else if (name == "seed") {
            if (!parse_integer(value, 0, 4294967295LL, seed)) {
                std::cerr << "error: --seed must be an integer from 0 to 4294967295\n";
                return 1;
            }
            have_seed = true;
        } else if (name == "end") {
            if (!parse_integer(value, 0, 1000000000, end) || !is_tick(static_cast<long>(end))) {
                std::cerr << "error: --end must be a checkpoint: " << CHECKPOINT_TEXT << "\n";
                return 1;
            }
        } else {
            std::cerr << "error: unknown option --" << name << "\n";
            usage();
            return 1;
        }
    }

    Job job;
    const bool any_setting = have_size || have_alpha || have_dens || have_seed;
    if (dir.empty() || end < 0 || (job_name.empty() && !any_setting)) {
        usage();
        return 1;
    }
    if (!job_name.empty()) {
        if (any_setting) {
            std::cerr << "error: give either --job or --size/--alpha/--dens/--seed, not both\n";
            return 1;
        }
        std::string error;
        if (!parse_job_id(job_name, job, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
    } else {
        if (!(have_size && have_alpha && have_dens && have_seed)) {
            std::cerr << "error: give all of --size, --alpha, --dens and --seed (or --job)\n";
            return 1;
        }
        job.size = static_cast<long>(size);
        job.alpha = alpha;
        job.dens = dens;
        job.seed = static_cast<uint32_t>(seed);
        job.id = job_id(job.size, job.alpha, job.dens, job.seed);
    }
    if (job.size < 4 || job.size * job.size > MAX_SITES) {
        std::cerr << "error: the lattice side must be from 4 to 10000\n";
        return 1;
    }
    if (!(job.alpha >= 0 && job.alpha <= 2)) {
        std::cerr << "error: alpha must be from 0 to 2\n";
        return 1;
    }
    if (!(job.dens >= 0) || !dens_supported(Model::Boson, job.dens)) {
        std::cerr << "error: dens must be from 0 to " << MAX_BOSON_DENS << "\n";
        return 1;
    }
    // The sweep's runners read the settings back from the folder's name, which allows only
    // digits and dots: an alpha or dens that JavaScript writes with an exponent (1e-7) has
    // no job name.
    Job named;
    std::string name_error;
    if (!parse_job_id(job.id, named, name_error)) {
        std::cerr << "error: the sweep's job names cannot hold these settings (the name would be "
                  << job.id << ")\n";
        return 1;
    }
    // The state files are little-endian; so is every machine this program targets.
    const uint32_t one = 1;
    if (*reinterpret_cast<const unsigned char*>(&one) != 1) {
        std::cerr << "error: the run mode writes little-endian lattices and needs a"
                     " little-endian machine\n";
        return 1;
    }

    try {
        return run_job(job, static_cast<long>(end), dir);
    } catch (const std::exception& e) {
        std::cerr << "error: " << job.id << ": " << e.what() << "\n";
        return 1;
    }
}

}  // namespace lgca
