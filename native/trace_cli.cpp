// native/trace_cli.cpp — the native program's `trace` mode: replay a sweep job from its
// seed and measure its island at every sample step.
//
// Usage:
//
//   lgca_native trace --job ID --end T [--check DIR] [--states DIR] OUT
//   lgca_native trace --size L --alpha A --dens D --seed S --end T [...] OUT
//
// The job is a sweep job (sweep_job.hpp; the settings of `lgca_native run`). The program
// runs it from its seed to step T with the same engine, so it passes through the same
// lattices as `run`. At each sample step (every 25 steps up to 5000, then every 100) it
// measures the lattice and adds a row to OUT, a JSON file:
//
//   {"job":{...}, "engineHash":..., "engine":{...}, "reached":T, "seconds":...,
//    "checked":{"dir":DIR,"rows":R,"snaps":S,"sha256":F} or null,
//    "columns":["t","N",...], "rows":[[0,19525,...],...],
//    "snaps":[{"t":0,"sha256":...,"file":...,"rng":[...]},...]}
//
// The columns are (measure.hpp): t; N, nmax, kmax, occupied, the cluster threshold and the
// heaviest cluster's mass and area (the sweeps' series fields, factor 2); then the island
// (Island: the busiest site and every occupied site connected to it): seed, islandMass,
// islandArea, islandPolar, wraps (1 if the island closes around the periodic box), ci, cj,
// rg, components (how many islands the lattice has) and largestMass; then, added on
// 8 October 2026 after those, largestArea (the heaviest island's sites), densestMass and
// densestArea (the island with the most particles per site), densest2Mass and densest2Area
// (the densest island of at least two sites). islandPolar is rounded to 6 decimals; ci, cj
// and rg to 3. A number that is undefined (the centre of an island that wraps, the densest
// island of two sites when every island is a single site) is written null. A complete
// trace has a row at every sample step from 0 to "reached", in order.
//
// `snaps` holds, at every checkpoint step (0 100 500 1000 2000 5000 10000 20000 40000
// 60000 80000, then every 20000 up to 2560000) up to T, the SHA-256 of the lattice, the
// name of its state file ("file", with --states) and the 33 generator words.
//
// --check DIR   DIR is a folder of the same job written by `lgca_native run` or by the
//               sweep's JavaScript runner. Its run.json must have the pinned engineHash
//               and the job's settings (id, size, alpha, dens, seed, sens, field; model
//               boson, alignment polar), as `run` requires to resume it, and hold a snap
//               at every checkpoint up to its last one and a series row at every sample
//               step from 0, at least to that checkpoint (a folder cut back to an earlier
//               checkpoint may keep later rows), each once and in order. Every series
//               row it saved up to T must be, character for character, the row this
//               replay measures, and every checkpoint it saved up to T must have the
//               replay's lattice SHA-256 and generator words. The first difference
//               stops the program with an error. "checked" names DIR, counts the rows
//               and checkpoints compared (all of DIR's up to T) and holds a fingerprint
//               ("sha256") of what was compared and of the trace it was compared with:
//               the SHA-256 of four lines, (1) the SHA-256 of DIR's compared rows and
//               snaps, their text as DIR saved it, each followed by a newline, rows
//               first; (2) the job as this program writes it ({"id":...,"field":"both"});
//               (3) OUT's rows and (4) OUT's snaps, each entry's text as written in OUT,
//               the entries joined by commas (as this program writes the list between
//               its brackets; space between entries is not part of it). DIR is only read.
// --states DIR  also save the lattice at each checkpoint, as DIR/state-t<T>.bin.gz (the
//               format of `run`), and name it in the snap ("file"). OUT may be in DIR,
//               but then its name must be plain ASCII and not run.json or a state
//               file's (or a state file's temporary name), in any case (a file system
//               may read two spellings as one name): that is refused before any file is
//               written (the folders may have been created).
//
// An existing OUT: it must be a complete trace of the same job (its settings, the
// engine, the columns, a row at every sample step and a snap at every checkpoint up to
// "reached", and "checked" well formed). Each row's values must be what this program
// writes: finite numbers, whole ones for the counts (t, N, nmax, kmax, occupied,
// threshold, clusterMass, clusterArea, seed from -1, islandMass, islandArea, wraps 0 or
// 1, components and the island masses and areas after it), null only for ci, cj and rg
// when wraps is 1 (and then all three), and for densest2Mass and densest2Area together
// (a number too large for a double, 1e999, reads as infinity and is refused). Anything
// else is refused with an error and left as it is (move it away to trace the job
// again). An OUT with other columns (written by
// an earlier version of this program) is refused in the same way: its rows cannot be
// continued. OUT is never replaced by a shorter trace: the program traces to the larger of
// T and OUT's "reached", and rewrites OUT only once the replay is past OUT's end, so a
// stopped program leaves OUT as it was. Where OUT and the replay overlap, every row and
// checkpoint of the replay must be OUT's, character for character (a difference stops
// the program). OUT counts as finished when it reached T, when (with --check) all of it
// counts as checked against DIR, and when (with --states) DIR holds the lattice of every
// one of its checkpoints and OUT names them; then the program says so and changes
// nothing. OUT counts as checked against DIR when its "checked" counts DIR's rows and
// checkpoints up to OUT's end and its fingerprint is that of DIR's rows and checkpoints
// and of OUT's own job, rows and snaps, as they are now (a change to either since the
// check shows). A trace written before the fingerprint was recorded, on 8 October 2026,
// has none, and the fingerprint's first form, the same day, covered only DIR's rows and
// checkpoints (line 1 alone): neither counts. Otherwise, with --check, the replay starts
// from the seed, so that every row is checked.
//
// Resume: with --states, the replay continues from the last checkpoint up to which DIR
// holds the lattice of every checkpoint (loaded and checked against OUT's SHA-256), or
// from the seed when there is none or OUT does not count as checked; without --states it
// starts from the seed. Missing lattices are saved again on the way.
//
// The states folder: before anything is saved there, the folder may hold no run.json (a
// run folder: its lattices are that run's) and no trace.json other than OUT, and its
// state files must all be this trace's checkpoints: those up to OUT's end must hold
// OUT's lattices; past it there may be only the next checkpoint's, which a stopped trace
// saved before it could write OUT (without an OUT: only a state-t0.bin.gz holding the
// job's starting lattice); that one is checked when the replay reaches it, and left as
// it is when the replay ends before it. A state file that is already there must hold the replay's
// lattice at that step (it is then kept as it is); one that holds another lattice is
// never overwritten: the program stops with an error, before it writes anything.
//
// One program at a time: the program locks OUT's folder and the states folder (flock on
// each folder, as `run` locks its folder) before it reads anything there, and holds the
// locks until it ends; a second trace or run that would work in either folder stops at
// once with an error. So one trace per folder: two traces whose OUTs share a folder do
// not run at the same time.
//
// OUT and the state files are replaced atomically (write_file_atomic), OUT at every
// checkpoint and at the end, so a stopped trace leaves complete files. A stop during a
// write can leave a stray temporary file; the next trace that goes on to write removes
// those of OUT and of the state files first (never the file OUT is or leads to, through
// a link). Prints "<ID>: t=<T>, <seconds>s" at the end.
// Errors stop the program with a message and exit status 1.

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
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

// Keep in step with COLUMNS in trace_sweep.py. New columns go at the end, so that the
// earlier ones keep their positions and their text.
const char* const COLUMNS[] = {
    "t",    "N",          "nmax",       "kmax",        "occupied", "threshold",
    "clusterMass", "clusterArea", "seed", "islandMass", "islandArea", "islandPolar",
    "wraps", "ci",        "cj",         "rg",          "components", "largestMass",
    "largestArea", "densestMass", "densestArea", "densest2Mass", "densest2Area"};
constexpr std::size_t COLUMN_COUNT = sizeof COLUMNS / sizeof COLUMNS[0];

// What each column holds, as row_json writes it: a whole number from 0 (seed from -1,
// wraps 0 or 1); islandPolar, any finite number; ci, cj and rg, a finite number, or null
// exactly when the island wraps; densest2Mass and densest2Area, whole numbers, or both null.
enum class Kind { Count, Seed, Flag, Number, Centre, Densest2 };
constexpr Kind KINDS[] = {
    Kind::Count, Kind::Count,  Kind::Count,  Kind::Count,  Kind::Count,  Kind::Count,
    Kind::Count, Kind::Count,  Kind::Seed,   Kind::Count,  Kind::Count,  Kind::Number,
    Kind::Flag,  Kind::Centre, Kind::Centre, Kind::Centre, Kind::Count,  Kind::Count,
    Kind::Count, Kind::Count,  Kind::Count,  Kind::Densest2, Kind::Densest2};
static_assert(sizeof KINDS / sizeof KINDS[0] == COLUMN_COUNT, "a kind for every column");
constexpr std::size_t WRAPS = 12, DENSEST2_MASS = 21, DENSEST2_AREA = 22;
constexpr double MAX_COUNT = std::numeric_limits<double>::max();   // any finite whole number

// Whether a row's values are what row_json writes; if not, the first column that is not
// (its name, in `bad`).
bool row_ok(const JsonValue& row, std::string& bad) {
    const auto is_null = [](const JsonValue& v) { return v.type == JsonValue::Type::Null; };
    const bool wraps = row.items[WRAPS].is_number() && row.items[WRAPS].number == 1;
    for (std::size_t c = 0; c < COLUMN_COUNT; ++c) {
        const JsonValue& v = row.items[c];
        bool ok = false;
        switch (KINDS[c]) {
            case Kind::Count: ok = is_integer(&v, 0, MAX_COUNT); break;
            case Kind::Seed: ok = is_integer(&v, -1, MAX_COUNT); break;
            case Kind::Flag: ok = is_integer(&v, 0, 1); break;
            case Kind::Number: ok = v.is_number() && std::isfinite(v.number); break;
            case Kind::Centre: ok = wraps ? is_null(v) : v.is_number() && std::isfinite(v.number); break;
            case Kind::Densest2:
                ok = is_null(row.items[DENSEST2_MASS]) && is_null(row.items[DENSEST2_AREA])
                         ? true
                         : is_integer(&v, 0, MAX_COUNT);
                break;
        }
        if (!ok) {
            bad = COLUMNS[c];
            return false;
        }
    }
    return true;
}

const std::string BUSY = ": another program (lgca_native run or trace) is working in this folder;"
                         " run this again when it has finished";

// v rounded to `digits` decimals (NaN stays NaN).
double rounded(double v, int digits) {
    const double scale = std::pow(10.0, digits);
    return std::round(v * scale) / scale;
}

// One row of OUT: the COLUMNS' values at step t.
std::string row_json(long t, const Measurements& m, const Island& s) {
    const double values[] = {
        static_cast<double>(t), m.N, static_cast<double>(m.nmax), static_cast<double>(m.kmax),
        static_cast<double>(m.occupied), m.threshold, m.clusterMass,
        static_cast<double>(m.clusterArea), static_cast<double>(s.seed), s.mass,
        static_cast<double>(s.area), rounded(s.polar, 6), s.wraps ? 1.0 : 0.0,
        rounded(s.ci, 3), rounded(s.cj, 3), rounded(s.rg, 3),
        static_cast<double>(s.components), s.largestMass,
        static_cast<double>(s.largestArea), s.densestMass, static_cast<double>(s.densestArea),
        s.densest2Mass, island_densest2_area(s)};
    static_assert(sizeof values / sizeof values[0] == COLUMN_COUNT);
    std::string out = "[";
    for (std::size_t c = 0; c < COLUMN_COUNT; ++c) {
        if (c) out += ",";
        out += json_number(values[c]);
    }
    return out + "]";
}

std::string snap_entry(long t, const std::string& sha, const uint32_t rng[33],
                       const std::string& file) {
    std::string s = "{\"t\":" + json_number(static_cast<double>(t)) + ",\"sha256\":" + json_string(sha);
    if (!file.empty()) s += ",\"file\":" + json_string(file);
    s += ",\"rng\":[";
    for (int q = 0; q < 33; ++q) {
        if (q) s += ",";
        s += std::to_string(rng[q]);
    }
    return s + "]}";
}

std::string job_json(const Job& job) {
    return "{\"id\":" + json_string(job.id) + ",\"size\":" + json_number(static_cast<double>(job.size)) +
           ",\"alpha\":" + json_number(job.alpha) + ",\"dens\":" + json_number(job.dens) +
           ",\"seed\":" + json_number(job.seed) + ",\"sens\":" + json_number(SENS) +
           ",\"field\":\"both\"}";
}

// The 33 generator words of a snap's "rng" (32 state words, then the index, below 32).
bool read_rng(const JsonValue* rng, uint32_t words[33]) {
    if (!rng || !rng->is_array() || rng->items.size() != 33) return false;
    for (int q = 0; q < 33; ++q) {
        if (!is_integer(&rng->items[q], 0, q < 32 ? 4294967295.0 : 31.0)) return false;
        words[q] = static_cast<uint32_t>(rng->items[q].number);
    }
    return true;
}

// How many of `steps` (in increasing order) are at most t.
std::size_t count_through(const std::vector<long>& steps, long t) {
    return static_cast<std::size_t>(std::upper_bound(steps.begin(), steps.end(), t) - steps.begin());
}

// ---------------------------------------------------------------------------------------
// The saved run a replay is checked against (--check DIR)
// ---------------------------------------------------------------------------------------

struct Saved {
    std::vector<long> row_t;                     // every sample step from 0, in order, at
                                                 // least to `last`
    std::vector<std::string> row_text;           // the rows' text as written
    struct Snap {
        long t;
        std::string sha;
        uint32_t rng[33];
        std::string text;
    };
    std::vector<Snap> snaps;   // every checkpoint to `last`, in order
    std::vector<long> snap_t;
    long last = -1;            // its last checkpoint

    std::size_t rows_through(long t) const { return count_through(row_t, t); }
    std::size_t snaps_through(long t) const { return count_through(snap_t, t); }

    // The fingerprint of the rows and snaps up to step t: the SHA-256 of their text as
    // saved, each followed by a newline, rows first.
    std::string fingerprint(long t) const {
        std::string all;
        for (std::size_t k = 0; k < rows_through(t); ++k) all += row_text[k] + "\n";
        for (std::size_t k = 0; k < snaps_through(t); ++k) all += snaps[k].text + "\n";
        return sha256_hex(all.data(), all.size());
    }
};

Saved load_saved(const std::string& dir, const Job& job) {
    namespace fs = std::filesystem;
    const std::string path = (fs::path(dir) / "run.json").string();
    std::vector<unsigned char> bytes;
    std::string error;
    // The limit `run` reads it with (run_cli.cpp): a run.json to step 2560000, the last
    // checkpoint, is about 12 MB.
    if (!read_file(path, std::size_t(64) << 20, bytes, error))
        throw std::runtime_error(path + ": cannot read: " + error);
    const std::string text(bytes.begin(), bytes.end());
    JsonValue run;
    if (!parse_json(text, run, error)) throw std::runtime_error(path + ": not valid JSON: " + error);
    const auto fail = [&path](const std::string& why) { return std::runtime_error(path + ": " + why); };
    if (!run.is_object()) throw fail("not a JSON object");
    const JsonValue* engine_hash = run.get("engineHash");
    if (!engine_hash || !engine_hash->is_string() || engine_hash->text != ENGINE_HASH)
        throw fail("engineHash is not the pinned engine's");
    const JsonValue* saved_job = run.get("job");
    const JsonValue* id = saved_job ? saved_job->get("id") : nullptr;
    if (!id || !id->is_string() || id->text != job.id) throw fail("not a run of " + job.id);
    // Its settings, as `run` checks them when it resumes a job.
    const auto same_number = [saved_job](const char* key, double value) {
        const JsonValue* v = saved_job->get(key);
        return v && v->is_number() && v->number == value;
    };
    const JsonValue* field = saved_job->get("field");
    if (!same_number("size", static_cast<double>(job.size)) || !same_number("alpha", job.alpha) ||
        !same_number("dens", job.dens) || !same_number("seed", job.seed) || !same_number("sens", SENS) ||
        !field || !field->is_string() || field->text != "both")
        throw fail("its job settings (size, alpha, dens, seed, sens, field) are not those of " + job.id);
    const JsonValue* model = run.get("model");
    const JsonValue* alignment = run.get("alignment");
    if (!model || !model->is_string() || model->text != "boson" || !alignment || !alignment->is_string() ||
        alignment->text != "polar")
        throw fail("its model and alignment are not boson and polar, those of " + job.id);

    Saved saved;
    // The snaps: every checkpoint from 0 to the last one, each once, in order.
    const JsonValue* snaps = run.get("snaps");
    if (!snaps || !snaps->is_array() || snaps->items.empty())
        throw fail("no checkpoint (\"snaps\" is missing or empty)");
    for (std::size_t k = 0; k < snaps->items.size(); ++k) {
        const JsonValue& snap = snaps->items[k];
        const JsonValue* t = snap.get("t");
        const JsonValue* sha = snap.get("sha256");
        const long expected = k < sizeof TICKS / sizeof TICKS[0] ? TICKS[k] : -1;
        if (!snap.is_object() || !is_integer(t, 0, 1e9) || t->number != static_cast<double>(expected))
            throw fail(std::string("its \"snaps\" must be every checkpoint step from 0 (") +
                       CHECKPOINT_TEXT + "), each once and in order: snap " +
                       std::to_string(k + 1) + " is not at t=" +
                       (expected < 0 ? std::string("a checkpoint") : std::to_string(expected)));
        Saved::Snap s;
        s.t = expected;
        if (!sha || !sha->is_string() || !is_lower_hex(sha->text, 64))
            throw fail("the snap at t=" + std::to_string(s.t) + " has no \"sha256\" of 64 lowercase hex digits");
        if (!read_rng(snap.get("rng"), s.rng))
            throw fail("the snap at t=" + std::to_string(s.t) + " has no valid \"rng\" (33 generator words)");
        s.sha = sha->text;
        s.text = text.substr(snap.begin, snap.end - snap.begin);
        saved.snaps.push_back(std::move(s));
        saved.snap_t.push_back(expected);
    }
    saved.last = saved.snap_t.back();

    // The series rows: every sample step from 0, each once, in order, at least to the last
    // checkpoint. A run cut back to an earlier checkpoint keeps its later rows (a resume
    // drops them): they are on the same schedule and are compared too.
    const JsonValue* series = run.get("series");
    if (!series || !series->is_array()) throw fail("no \"series\" array");
    long expected = 0;
    for (std::size_t k = 0; k < series->items.size(); ++k) {
        const JsonValue& row = series->items[k];
        const JsonValue* t = row.get("t");
        if (!row.is_object() || !is_integer(t, 0, 1e9) || t->number != static_cast<double>(expected))
            throw fail("its \"series\" rows must be every sample step from 0 (every 25 steps to"
                       " 5000, then every 100), each once and in order, at least to its last"
                       " checkpoint, t=" + std::to_string(saved.last) + ": row " + std::to_string(k + 1) +
                       " is not at t=" + std::to_string(expected));
        saved.row_t.push_back(expected);
        saved.row_text.push_back(text.substr(row.begin, row.end - row.begin));
        expected += expected < 5000 ? 25 : 100;
    }
    if (saved.row_t.empty() || saved.row_t.back() < saved.last)
        throw fail("its \"series\" rows must be every sample step from 0, at least to its last"
                   " checkpoint, t=" + std::to_string(saved.last) + ": they stop before t=" +
                   std::to_string(expected));
    return saved;
}

// ---------------------------------------------------------------------------------------
// An existing OUT
// ---------------------------------------------------------------------------------------

struct Existing {
    bool exists = false;
    long reached = -1;
    double seconds = 0;
    std::vector<std::string> row_text;          // every sample step to `reached`, as written
    struct Snap {
        long t;
        std::string sha;
        uint32_t rng[33];
        bool file;          // names its state file
        std::string text;   // as written
    };
    std::vector<Snap> snaps;   // every checkpoint to `reached`
    bool checked = false;
    std::string checked_dir, checked_sha;   // checked_sha "": written before it was recorded
    long checked_rows = 0, checked_snaps = 0;
};

// Read OUT if it exists. Throws std::runtime_error, with the reason, if it is not a
// complete trace of this job with the current columns.
Existing read_existing(const std::string& out, const Job& job) {
    namespace fs = std::filesystem;
    Existing ex;
    std::error_code ec;
    const fs::file_status status = fs::symlink_status(out, ec);
    if (status.type() == fs::file_type::not_found) return ex;
    if (ec) throw std::runtime_error(out + ": " + ec.message());
    ex.exists = true;
    std::vector<unsigned char> bytes;
    std::string error;
    // A trace to step 2560000, the last checkpoint, holds 25751 rows and 135 snaps, about
    // 4 MB at any lattice size; the limit keeps a damaged file from filling memory.
    if (!read_file(out, std::size_t(64) << 20, bytes, error))
        throw std::runtime_error(out + ": cannot read: " + error);
    const std::string text(bytes.begin(), bytes.end());
    JsonValue old;
    if (!parse_json(text, old, error)) throw std::runtime_error(out + ": not valid JSON: " + error);
    const JsonValue* old_job = old.get("job");
    const JsonValue* id = old_job ? old_job->get("id") : nullptr;
    const JsonValue* engine_hash = old.get("engineHash");
    if (!old.is_object() || !id || !id->is_string() || id->text != job.id || !engine_hash ||
        !engine_hash->is_string() || engine_hash->text != ENGINE_HASH)
        throw std::runtime_error(out + ": exists and is not a trace of " + job.id +
                                 "; refusing to replace it");
    // Other columns (an earlier version of this program): neither continued (its rows would
    // mix two layouts) nor replaced (a stopped replacement would lose it).
    const JsonValue* columns = old.get("columns");
    bool other_columns = !columns || !columns->is_array() || columns->items.size() != COLUMN_COUNT;
    for (std::size_t c = 0; !other_columns && c < COLUMN_COUNT; ++c)
        other_columns = !columns->items[c].is_string() || columns->items[c].text != COLUMNS[c];
    if (other_columns)
        throw std::runtime_error(out + " has other columns (written by an earlier version);"
                                 " it is not continued or replaced: write the trace elsewhere");
    const auto invalid = [&out, &job](const std::string& why) {
        return std::runtime_error(out + ": exists but is not a complete trace of " + job.id + " (" + why +
                                  "); refusing to replace it: move it away to trace the job again");
    };

    // The job's settings, as its name gives them.
    const auto same_number = [old_job](const char* key, double value) {
        const JsonValue* v = old_job->get(key);
        return v && v->is_number() && v->number == value;
    };
    const JsonValue* field = old_job->get("field");
    if (!same_number("size", static_cast<double>(job.size)) || !same_number("alpha", job.alpha) ||
        !same_number("dens", job.dens) || !same_number("seed", job.seed) || !same_number("sens", SENS) ||
        !field || !field->is_string() || field->text != "both")
        throw invalid("its job settings are not those of " + job.id);

    const JsonValue* reached = old.get("reached");
    if (!is_integer(reached, 0, 1e9) || !is_tick(static_cast<long>(reached->number)))
        throw invalid("\"reached\" is not a checkpoint step");
    ex.reached = static_cast<long>(reached->number);
    const JsonValue* seconds = old.get("seconds");
    if (!seconds || !seconds->is_number() || !std::isfinite(seconds->number) || seconds->number < 0)
        throw invalid("\"seconds\" is not a number of seconds");
    ex.seconds = seconds->number;

    // A row at every sample step up to `reached`, each of COLUMN_COUNT values as row_json
    // writes them.
    const JsonValue* rows = old.get("rows");
    if (!rows || !rows->is_array()) throw invalid("no \"rows\" array");
    const std::vector<long> samples = sample_steps(ex.reached);
    if (rows->items.size() != samples.size())
        throw invalid(std::to_string(rows->items.size()) + " rows where every sample step to t=" +
                      std::to_string(ex.reached) + " makes " + std::to_string(samples.size()));
    for (std::size_t k = 0; k < samples.size(); ++k) {
        const JsonValue& row = rows->items[k];
        if (!row.is_array() || row.items.size() != COLUMN_COUNT || !row.items[0].is_number() ||
            row.items[0].number != static_cast<double>(samples[k]))
            throw invalid("row " + std::to_string(k + 1) + " is not " + std::to_string(COLUMN_COUNT) +
                          " values starting with its step, t=" + std::to_string(samples[k]));
        std::string bad;
        if (!row_ok(row, bad))
            throw invalid("row " + std::to_string(k + 1) + ", t=" + std::to_string(samples[k]) +
                          ", has an invalid \"" + bad + "\": each value is a finite number, a whole one"
                          " where it counts something; null only for ci, cj and rg of an island that"
                          " wraps, and for densest2Mass and densest2Area together");
        ex.row_text.push_back(text.substr(row.begin, row.end - row.begin));
    }

    // A snap at every checkpoint up to `reached`.
    const JsonValue* snaps = old.get("snaps");
    if (!snaps || !snaps->is_array()) throw invalid("no \"snaps\" array");
    const std::vector<long> ticks = tick_steps(ex.reached);
    if (snaps->items.size() != ticks.size())
        throw invalid(std::to_string(snaps->items.size()) + " snaps where every checkpoint to t=" +
                      std::to_string(ex.reached) + " makes " + std::to_string(ticks.size()));
    for (std::size_t k = 0; k < ticks.size(); ++k) {
        const JsonValue& snap = snaps->items[k];
        const JsonValue* t = snap.get("t");
        const JsonValue* sha = snap.get("sha256");
        const JsonValue* file = snap.get("file");
        Existing::Snap s;
        s.t = ticks[k];
        s.file = file != nullptr;
        if (!snap.is_object() || !t || !t->is_number() || t->number != static_cast<double>(s.t) || !sha ||
            !sha->is_string() || !is_lower_hex(sha->text, 64) || !read_rng(snap.get("rng"), s.rng) ||
            (file && (!file->is_string() || file->text != state_file_name(s.t))))
            throw invalid("snap " + std::to_string(k + 1) + " is not the checkpoint t=" +
                          std::to_string(s.t) + " with its \"sha256\", 33 \"rng\" words and"
                          " (if any) \"file\" " + state_file_name(s.t));
        s.sha = sha->text;
        s.text = text.substr(snap.begin, snap.end - snap.begin);
        ex.snaps.push_back(s);
    }

    // "checked": null, or what was compared with a saved run.
    const JsonValue* checked = old.get("checked");
    if (!checked) throw invalid("no \"checked\"");
    if (checked->is_object()) {
        ex.checked = true;
        const JsonValue* dir = checked->get("dir");
        const JsonValue* cr = checked->get("rows");
        const JsonValue* cs = checked->get("snaps");
        const JsonValue* sha = checked->get("sha256");
        if (!dir || !dir->is_string() || !is_integer(cr, 0, static_cast<double>(samples.size())) ||
            !is_integer(cs, 0, static_cast<double>(ticks.size())) ||
            (sha && (!sha->is_string() || !is_lower_hex(sha->text, 64))))
            throw invalid("\"checked\" is not {\"dir\", \"rows\", \"snaps\"[, \"sha256\"]} within its rows"
                          " and snaps");
        ex.checked_dir = dir->text;
        ex.checked_rows = static_cast<long>(cr->number);
        ex.checked_snaps = static_cast<long>(cs->number);
        if (sha) ex.checked_sha = sha->text;
    } else if (checked->type != JsonValue::Type::Null) {
        throw invalid("\"checked\" is neither null nor an object");
    }
    return ex;
}

// The fingerprint "checked" records (the top of this file): the SHA-256 of four lines, the
// saved run's fingerprint up to the trace's end, the job, and the trace's rows and snaps
// (each entry's text, the entries joined by commas).
std::string check_fingerprint(const std::string& saved_fingerprint, const Job& job, const std::string& rows,
                              const std::string& snaps) {
    const std::string all = saved_fingerprint + "\n" + job_json(job) + "\n" + rows + "\n" + snaps + "\n";
    return sha256_hex(all.data(), all.size());
}

// Whether all of OUT counts as checked against the saved run in `dir` (the top of this
// file); if not, why, in `why`.
bool covers(const Existing& ex, const Saved& saved, const Job& job, const std::string& dir, std::string& why) {
    if (!ex.checked) {
        why = "it was not checked";
        return false;
    }
    const long T = ex.reached;
    const std::size_t rows = saved.rows_through(T), snaps = saved.snaps_through(T);
    if (static_cast<std::size_t>(ex.checked_rows) != rows ||
        static_cast<std::size_t>(ex.checked_snaps) != snaps) {
        why = "it was checked against " + std::to_string(ex.checked_rows) + " rows and " +
              std::to_string(ex.checked_snaps) + " checkpoints, and " + dir + " has " + std::to_string(rows) +
              " and " + std::to_string(snaps) + " up to t=" + std::to_string(T);
        return false;
    }
    if (ex.checked_sha.empty()) {
        why = "it was checked before checks were fingerprinted, so nothing shows that " + dir +
              " still holds what it was checked against";
        return false;
    }
    std::string own_rows, own_snaps;
    for (const std::string& row : ex.row_text) append_item(own_rows, row);
    for (const Existing::Snap& s : ex.snaps) append_item(own_snaps, s.text);
    const std::string source = saved.fingerprint(T);
    if (ex.checked_sha == check_fingerprint(source, job, own_rows, own_snaps)) return true;
    if (ex.checked_sha == source)
        why = "its check was recorded before it covered the trace's own rows and checkpoints, so nothing"
              " shows that they are the ones that were checked";
    else
        why = "its rows and checkpoints, or those of " + dir + " up to t=" + std::to_string(T) +
              ", are not the ones that were checked";
    return false;
}

// ---------------------------------------------------------------------------------------
// The states folder
// ---------------------------------------------------------------------------------------

// Whether `path` holds the lattice whose SHA-256 is `sha`: true if so, false if there is
// no file; throws if it holds another lattice or cannot be read (it is never overwritten).
bool holds_lattice(const std::string& path, std::size_t lattice_bytes, const std::string& sha, long t) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (fs::symlink_status(path, ec).type() == fs::file_type::not_found) return false;
    std::vector<unsigned char> bytes;
    std::string error;
    if (!read_file(path, lattice_bytes, bytes, error)) {
        throw std::runtime_error(path + ": cannot be read (" + error + "); it is not overwritten: move it"
                                 " away to save this trace's lattice at t=" + std::to_string(t));
    }
    if (bytes.size() != lattice_bytes || sha256_hex(bytes.data(), bytes.size()) != sha)
        throw std::runtime_error(path + ": holds another lattice than this trace's at t=" + std::to_string(t) +
                                 "; it is not overwritten: move it away, or give another --states folder");
    return true;
}

// Before anything is saved in `states` (the top of this file): no run.json, no trace.json
// but OUT, and only state files of this trace's checkpoints: those up to OUT's end (each
// checked against OUT's SHA-256 by states_prefix), and past it only the one a stopped
// trace saved before it could write OUT: one file, at the checkpoint after OUT's end,
// which the replay checks before it writes there (and leaves alone if it stops short of
// it); without an OUT, a lone state-t0.bin.gz that holds this job's starting lattice
// (`start_sha`).
void check_state_folder(const std::string& states, const std::string& out, const Existing& ex,
                        const std::string& start_sha, std::size_t lattice_bytes) {
    namespace fs = std::filesystem;
    std::error_code ec;
    const fs::path run = fs::path(states) / "run.json", trace = fs::path(states) / "trace.json";
    if (fs::symlink_status(run, ec).type() != fs::file_type::not_found)
        throw std::runtime_error(states + ": holds a run.json: it is a run's folder (lgca_native run),"
                                 " and its lattices are that run's; save the trace's lattices in another"
                                 " folder");
    if (fs::symlink_status(trace, ec).type() != fs::file_type::not_found &&
        !(ex.exists && fs::equivalent(trace, out, ec)))
        throw std::runtime_error(states + ": holds another trace (trace.json, not " + out +
                                 "), whose lattices these are; save this trace's lattices in another"
                                 " folder");
    long next = -1;   // the checkpoint after OUT's end (0 without an OUT)
    for (const long tick : TICKS)
        if (tick > ex.reached) {
            next = tick;
            break;
        }
    for (const std::string& name : folder_entries(states)) {
        long t = -1;
        if (!is_state_file_name(name, &t)) continue;
        if (t < 0 || !is_tick(t) || state_file_name(t) != name)
            throw std::runtime_error(
                states + ": holds " + name + ", which is not a checkpoint's state file: refusing to mix"
                " another program's files with this trace's lattices; give another --states folder, or"
                " move the file away");
        if (t <= ex.reached || t == next) continue;
        throw std::runtime_error(
            states + ": holds " + name + ", but " +
            (ex.exists ? out + " ends at t=" + std::to_string(ex.reached) : out + " does not exist") +
            ", so no trace names it: refusing to mix it with this trace's lattices; move it away (or put"
            " back the trace that names it)");
    }
    if (!ex.exists) {
        const std::string path = (fs::path(states) / state_file_name(0)).string();
        holds_lattice(path, lattice_bytes, start_sha, 0);   // throws if it is another lattice
    }
}

// The last checkpoint up to which `states` holds OUT's lattice at every checkpoint (-1:
// not even step 0's). Every state file that is there, after a gap too, must hold OUT's
// lattice (holds_lattice).
long states_prefix(const Existing& ex, const std::string& states, std::size_t lattice_bytes) {
    namespace fs = std::filesystem;
    long prefix = -1;
    bool gap = false;
    for (const Existing::Snap& s : ex.snaps) {
        const bool there = holds_lattice((fs::path(states) / state_file_name(s.t)).string(), lattice_bytes,
                                         s.sha, s.t);
        if (!there) gap = true;
        else if (!gap) prefix = s.t;
    }
    return prefix;
}

// Remove the stray temporary files (atomic_temp_target) in `dir` of the files `ours`
// accepts, but never the file OUT is or leads to (a link can make OUT a file of that
// name). The caller holds the folder's lock, so no other program is writing them.
void remove_stray_temp_files(const std::string& dir, const std::string& out,
                             const std::function<bool(const std::string&)>& ours) {
    namespace fs = std::filesystem;
    for (const std::string& name : folder_entries(dir)) {
        const std::string target = atomic_temp_target(name);
        if (target.empty() || !ours(target)) continue;
        const fs::path path = fs::path(dir) / name;
        std::error_code ec;
        // A folder, or a link to one (OUT's path may pass through it), is not a stray file.
        if (fs::is_directory(fs::status(path, ec)) || fs::equivalent(path, out, ec)) continue;
        if (!fs::remove(path, ec) && ec)
            throw std::runtime_error(path.string() + ": cannot remove this stray temporary file: " +
                                     ec.message());
    }
}

// ---------------------------------------------------------------------------------------

std::string trace_json(const Job& job, long reached, double seconds, const std::string& check_dir,
                       long checked_rows, long checked_snaps, const std::string& fingerprint,
                       const std::string& rows, const std::string& snaps) {
    std::string s;
    s.reserve(rows.size() + snaps.size() + 2048);
    s += "{\"job\":" + job_json(job);
    s += ",\"engineHash\":" + json_string(ENGINE_HASH);
    s += ",\"engine\":" + engine_json("lgca_native trace");
    s += ",\"reached\":" + json_number(static_cast<double>(reached));
    s += ",\"seconds\":" + json_number(seconds);
    s += ",\"checked\":";
    s += !check_dir.empty() ? "{\"dir\":" + json_string(check_dir) + ",\"rows\":" +
                                  std::to_string(checked_rows) + ",\"snaps\":" +
                                  std::to_string(checked_snaps) + ",\"sha256\":" +
                                  json_string(fingerprint) + "}"
                            : std::string("null");
    s += ",\"columns\":[";
    for (std::size_t c = 0; c < COLUMN_COUNT; ++c) {
        if (c) s += ",";
        s += json_string(COLUMNS[c]);
    }
    s += "],\"rows\":[" + rows + "],\"snaps\":[" + snaps + "]}";
    return s;
}

int trace_job(const Job& job, long end, const std::string& check_dir, const std::string& states,
              const std::string& out) {
    namespace fs = std::filesystem;
    std::error_code ec;
    // The folders, each locked before anything in it is read: OUT's, and the states folder.
    const fs::path out_path(out);
    const std::string out_name = out_path.filename().string();
    if (out_name.empty() || out_name == "." || out_name == "..")
        throw std::runtime_error(out + ": OUT must name a file");
    const std::string out_dir = out_path.has_parent_path() ? out_path.parent_path().string() : ".";
    const auto make_folder = [&ec](const std::string& dir) {
        fs::create_directories(dir, ec);
        if (ec || !fs::is_directory(dir))
            throw std::runtime_error(dir + ": cannot create the folder" + (ec ? ": " + ec.message() : ""));
    };
    make_folder(out_dir);
    bool same_folder = false;
    if (!states.empty()) {
        make_folder(states);
        same_folder = fs::equivalent(states, out_dir, ec);
        // In the states folder, OUT may not take a name that folder keeps for other files.
        // A file system can read two spellings as one name (macOS ignores case, and folds
        // some other letters into ASCII ones: a long s is an s), so there OUT's name must
        // be plain ASCII, and it is compared in lower case.
        std::string name = out_name;
        bool ascii = true;
        for (char& c : name) {
            ascii = ascii && static_cast<unsigned char>(c) < 0x80;
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        if (same_folder && (!ascii || name == "run.json" || is_state_file_name(name) ||
                            is_state_file_name(atomic_temp_target(name))))
            throw std::runtime_error(out + ": OUT is in the states folder " + states + ", where its name must"
                                     " be plain ASCII and not run.json or a state file's (state-t<T>.bin.gz,"
                                     " or its temporary name), in any case; name it otherwise, as trace.json");
    }
    const FolderLock out_lock(out_dir, out_dir + BUSY);
    std::optional<FolderLock> states_lock;
    if (!states.empty() && !same_folder) states_lock.emplace(states, states + BUSY);

    const Params params = job_params(job);
    std::optional<Lgca> engine;
    engine.emplace(params);   // seeded and filled: step 0
    Measurer measurer(params.W, params.H);
    const std::size_t cells = engine->cell_count();
    const std::size_t lattice_bytes = cells * sizeof(int32_t);

    const Existing ex = read_existing(out, job);
    // OUT is never replaced by a shorter trace.
    const long target = ex.exists ? std::max(end, ex.reached) : end;
    std::optional<Saved> saved;
    if (!check_dir.empty()) saved = load_saved(check_dir, job);
    if (!states.empty())
        check_state_folder(states, out, ex, sha256_hex(engine->occ(), lattice_bytes), lattice_bytes);

    std::string why_unchecked;
    const bool covered = !saved || (ex.exists && covers(ex, *saved, job, check_dir, why_unchecked));
    long prefix = -1;   // the states folder holds OUT's lattices of every checkpoint to here
    bool named = true;  // OUT's snaps name their state files
    if (!states.empty() && ex.exists) {
        prefix = states_prefix(ex, states, lattice_bytes);
        for (const Existing::Snap& s : ex.snaps) named = named && s.file;
    }
    if (ex.exists && ex.reached >= end && covered &&
        (states.empty() || (prefix == ex.reached && named))) {
        std::cout << job.id << ": already at t=" << ex.reached << ", nothing to do\n";
        return 0;
    }
    if (ex.exists && !covered)
        std::cout << job.id << ": " << out << " does not count as checked against " << check_dir << " ("
                  << why_unchecked << "): the replay starts from the seed\n";

    // Where the replay starts: the last checkpoint up to which every lattice is saved, when
    // OUT's rows count as checked; otherwise the seed.
    const long start = ex.exists && covered && !states.empty() ? prefix : -1;
    Lgca& sim = *engine;
    std::string rows, snaps, error;
    std::size_t row_index = 0, snap_index = 0;   // the replay's place in OUT's rows and snaps
    double previous_seconds = 0;
    long checked_rows = 0, checked_snaps = 0;
    if (start >= 0) {
        std::vector<unsigned char> bytes;
        const std::string path = (fs::path(states) / state_file_name(start)).string();
        const std::size_t k = count_through(tick_steps(start), start) - 1;
        if (!read_file(path, lattice_bytes, bytes, error) || bytes.size() != lattice_bytes ||
            sha256_hex(bytes.data(), bytes.size()) != ex.snaps[k].sha)
            throw std::runtime_error(path + ": cannot be read again, or has changed");
        std::memcpy(sim.occ(), bytes.data(), lattice_bytes);
        sim.set_rng_state(ex.snaps[k].rng);
        sim.set_t(static_cast<std::int64_t>(start));
        row_index = sample_steps(start).size();
        snap_index = k + 1;
        for (std::size_t r = 0; r < row_index; ++r) append_item(rows, ex.row_text[r]);
        for (std::size_t s = 0; s < snap_index; ++s)
            append_item(snaps, snap_entry(ex.snaps[s].t, ex.snaps[s].sha, ex.snaps[s].rng,
                                          state_file_name(ex.snaps[s].t)));
        previous_seconds = ex.seconds;
        if (saved) {
            checked_rows = static_cast<long>(saved->rows_through(start));
            checked_snaps = static_cast<long>(saved->snaps_through(start));
        }
    }

    // About to write: first clear what a stopped trace left.
    remove_stray_temp_files(out_dir, out, [&out_name](const std::string& target_name) {
        return target_name == out_name;
    });
    if (!states.empty())
        remove_stray_temp_files(states, out, [](const std::string& target_name) {
            return is_state_file_name(target_name);
        });

    const auto begin = std::chrono::steady_clock::now();
    const auto seconds = [&begin, previous_seconds]() {
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - begin).count();
        return previous_seconds + static_cast<double>(ms) / 1000;
    };
    long written = -1;
    const auto write_out = [&](long reached) {
        std::string fingerprint;
        if (saved) {
            // Every row and checkpoint the saved run has up to `reached` was compared.
            if (static_cast<std::size_t>(checked_rows) != saved->rows_through(reached) ||
                static_cast<std::size_t>(checked_snaps) != saved->snaps_through(reached))
                throw std::runtime_error("the replay compared " + std::to_string(checked_rows) + " rows and " +
                                         std::to_string(checked_snaps) + " checkpoints up to t=" +
                                         std::to_string(reached) + ", not all of " + check_dir + "'s");
            fingerprint = check_fingerprint(saved->fingerprint(reached), job, rows, snaps);
        }
        const std::string text = trace_json(job, reached, seconds(), check_dir, checked_rows,
                                            checked_snaps, fingerprint, rows, snaps);
        if (!write_file_atomic(out, text.data(), text.size(), error)) throw std::runtime_error(error);
        written = reached;
    };
    const auto differs = [&out](long t, const char* what) {
        return std::runtime_error("the replay differs from " + out + " at the " + what + " t=" +
                                  std::to_string(t) + ": it was not written by this job's replay with"
                                  " this engine; it is left as it is");
    };

    std::vector<unsigned char> gz;
    for (long t = start + 1; t <= target; ++t) {
        if (t) sim.step(SENS);
        if (!is_sample(t)) continue;
        const int32_t* occ = sim.occ();
        const Measurements m = measurer.measure(occ, 2);
        bool negative = false;
        for (std::size_t c = 0; c < cells; ++c) negative |= occ[c] < 0;
        if (m.N != static_cast<double>(sim.nuparts()) || negative)
            throw std::runtime_error("Invalid particle count at t=" + std::to_string(t));
        if (saved) {
            const std::size_t k = count_through(saved->row_t, t);
            if (k > 0 && saved->row_t[k - 1] == t) {
                if (saved->row_text[k - 1] != series_json(t, m))
                    throw std::runtime_error("the replay differs from " + check_dir +
                                             " at the series row t=" + std::to_string(t));
                ++checked_rows;
            }
        }
        const std::string row = row_json(t, m, measurer.island(occ));
        if (ex.exists && t <= ex.reached && row != ex.row_text[row_index]) throw differs(t, "row");
        ++row_index;
        append_item(rows, row);
        if (!is_tick(t)) continue;
        const std::string sha = sha256_hex(occ, lattice_bytes);
        uint32_t rng[33];
        sim.get_rng_state(rng);
        if (saved) {
            const std::size_t k = count_through(saved->snap_t, t);
            if (k > 0 && saved->snap_t[k - 1] == t) {
                if (saved->snaps[k - 1].sha != sha ||
                    std::memcmp(saved->snaps[k - 1].rng, rng, sizeof rng) != 0)
                    throw std::runtime_error("the replay differs from " + check_dir +
                                             " at the checkpoint t=" + std::to_string(t));
                ++checked_snaps;
            }
        }
        if (ex.exists && t <= ex.reached &&
            (ex.snaps[snap_index].sha != sha || std::memcmp(ex.snaps[snap_index].rng, rng, sizeof rng) != 0))
            throw differs(t, "checkpoint");
        ++snap_index;
        std::string file;
        if (!states.empty()) {
            file = state_file_name(t);
            const std::string path = (fs::path(states) / file).string();
            // A file already there must be this lattice (then it stays as it is).
            if (!holds_lattice(path, lattice_bytes, sha, t) &&
                (!gzip_compress(occ, lattice_bytes, gz, error) ||
                 !write_file_atomic(path, gz.data(), gz.size(), error)))
                throw std::runtime_error(error);
        }
        append_item(snaps, snap_entry(t, sha, rng, file));
        if (!ex.exists || t >= ex.reached) write_out(t);
    }
    // Nothing left to replay (OUT only lacked the names of its state files): rewrite it.
    if (written != target) write_out(target);
    char line[64];
    std::snprintf(line, sizeof line, "%.1f", seconds());
    std::cout << job.id << ": t=" << target << ", " << line << "s";
    if (saved) std::cout << " (checked " << checked_rows << " rows, " << checked_snaps << " checkpoints)";
    std::cout << "\n";
    return 0;
}

void usage() {
    std::cerr << "usage: lgca_native trace --job ID --end T [--check DIR] [--states DIR] OUT\n"
                 "       lgca_native trace --size L --alpha A --dens D --seed S --end T"
                 " [--check DIR] [--states DIR] OUT\n";
}

}  // namespace

int trace_main(int argc, char** argv) {
    std::string job_name, out, check_dir, states;
    bool have_size = false, have_alpha = false, have_dens = false, have_seed = false;
    long long size = 0, seed = 0, end = -1;
    double alpha = 0, dens = 0;

    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) {
            if (!out.empty()) {
                std::cerr << "error: give one output file (also got " << arg << ")\n";
                usage();
                return 1;
            }
            out = arg;
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
        } else if (name == "check") {
            check_dir = value;
        } else if (name == "states") {
            states = value;
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
    if (out.empty() || end < 0 || (job_name.empty() && !any_setting)) {
        usage();
        return 1;
    }
    std::string error;
    if (!job_name.empty()) {
        if (any_setting) {
            std::cerr << "error: give either --job or --size/--alpha/--dens/--seed, not both\n";
            return 1;
        }
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
        Job named;
        if (!parse_job_id(job.id, named, error)) {
            std::cerr << "error: the sweep's job names cannot hold these settings (the name would"
                         " be " << job.id << ")\n";
            return 1;
        }
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
    const uint32_t one = 1;
    if (*reinterpret_cast<const unsigned char*>(&one) != 1) {
        std::cerr << "error: the trace mode needs a little-endian machine\n";
        return 1;
    }

    try {
        return trace_job(job, static_cast<long>(end), check_dir, states, out);
    } catch (const std::exception& e) {
        std::cerr << "error: " << job.id << ": " << e.what() << "\n";
        return 1;
    }
}

}  // namespace lgca
