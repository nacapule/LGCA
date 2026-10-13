// native/sweep_job.hpp — what the `run` and `trace` modes share: a sweep job's fixed
// settings, its name, its sample and checkpoint steps, and small helpers for run.json.
//
// A sweep job is one run of the density and alpha sweep's JavaScript runner: boson
// model, power kernel with exponent alpha, field from the site and its 6 neighbours
// ("both"), polar alignment, sensitivity 6, an L x L lattice. Only the side L, alpha,
// the density and the seed vary; the job's name holds all four.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "json.hpp"
#include "lgca_engine.hpp"

namespace lgca {

// The SHA-256 of the JavaScript engine source that the density and alpha sweep's runner
// pinned and wrote as "engineHash". This program reproduces that engine's dynamics, so
// every native build writes the same value: it names the dynamics, not the build (the
// "engine" field, engine_json, names the build).
inline constexpr const char* ENGINE_HASH =
    "34edbb46d96ee9aaea500e8c5e0927a5858366a2ed52fc3db1284d37b7d5ff6d";

// The sweep's fixed settings.
inline constexpr double SENS = 6;
// The checkpoints: the steps at which a job saves its lattice, and at which it may end and
// later continue. Up to 80000 they are the JavaScript runner's; after it they go on every
// 20000 steps up to 2560000. Keep in step with CHECKPOINTS in sweep.py.
inline constexpr long TICKS[] = {0,       100,     500,     1000,    2000,    5000,    10000,
                                 20000,   40000,   60000,   80000,   100000,  120000,  140000,
                                 160000,  180000,  200000,  220000,  240000,  260000,  280000,
                                 300000,  320000,  340000,  360000,  380000,  400000,  420000,
                                 440000,  460000,  480000,  500000,  520000,  540000,  560000,
                                 580000,  600000,  620000,  640000,  660000,  680000,  700000,
                                 720000,  740000,  760000,  780000,  800000,  820000,  840000,
                                 860000,  880000,  900000,  920000,  940000,  960000,  980000,
                                 1000000, 1020000, 1040000, 1060000, 1080000, 1100000, 1120000,
                                 1140000, 1160000, 1180000, 1200000, 1220000, 1240000, 1260000,
                                 1280000, 1300000, 1320000, 1340000, 1360000, 1380000, 1400000,
                                 1420000, 1440000, 1460000, 1480000, 1500000, 1520000, 1540000,
                                 1560000, 1580000, 1600000, 1620000, 1640000, 1660000, 1680000,
                                 1700000, 1720000, 1740000, 1760000, 1780000, 1800000, 1820000,
                                 1840000, 1860000, 1880000, 1900000, 1920000, 1940000, 1960000,
                                 1980000, 2000000, 2020000, 2040000, 2060000, 2080000, 2100000,
                                 2120000, 2140000, 2160000, 2180000, 2200000, 2220000, 2240000,
                                 2260000, 2280000, 2300000, 2320000, 2340000, 2360000, 2380000,
                                 2400000, 2420000, 2440000, 2460000, 2480000, 2500000, 2520000,
                                 2540000, 2560000};
// The checkpoints in words, for messages.
inline constexpr const char* CHECKPOINT_TEXT =
    "0 100 500 1000 2000 5000 10000, then every 20000 up to 2560000";
static_assert(
    [] {
        constexpr std::size_t n = sizeof TICKS / sizeof TICKS[0];
        for (std::size_t k = 1; k < n; ++k)
            if (TICKS[k] <= TICKS[k - 1] || (TICKS[k] > 20000 && TICKS[k] - TICKS[k - 1] != 20000))
                return false;
        return TICKS[0] == 0 && TICKS[6] == 10000 && TICKS[7] == 20000 && TICKS[n - 1] == 2560000;
    }(),
    "TICKS must be the steps CHECKPOINT_TEXT names");

inline bool is_tick(long t) {
    for (const long tick : TICKS)
        if (t == tick) return true;
    return false;
}

// The runner's isSample = t <= 5000 ? t % 25 === 0 : t % 100 === 0.
inline bool is_sample(long t) { return t <= 5000 ? t % 25 == 0 : t % 100 == 0; }

// The sample steps from 0 to t, in order: what a complete `series` (run.json) or `rows`
// (trace.json) holds when the job has reached t.
std::vector<long> sample_steps(long t);

// The checkpoint steps from 0 to t, in order: a complete `snaps` at t.
std::vector<long> tick_steps(long t);

// A checkpoint's state file name, state-t<T>.bin.gz.
std::string state_file_name(long t);

// Whether `name` is a state file's name, state-t<digits>.bin.gz; if so and `t` is given,
// its step (or -1 if the digits do not fit in a long).
bool is_state_file_name(const std::string& name, long* t = nullptr);

inline constexpr long MAX_SITES = 100000000;   // the largest lattice: 10^8 sites (L = 10000)

// The job: the runner's `job` object {id, size, alpha, dens, seed, sens, field}.
struct Job {
    std::string id;
    long size = 0;
    double alpha = 0, dens = 0;
    uint32_t seed = 0;
};

// The engine settings of a job (seed and density included).
Params job_params(const Job& job);

// The whole text must be a decimal integer from min_value to max_value.
bool parse_integer(const std::string& text, long long min_value, long long max_value,
                   long long& value);

// The whole text must be a finite number.
bool parse_number(const std::string& text, double& value);

// The job name the sweep uses: both-L<size>-a<alpha>-d<dens>-seed<seed>.
std::string job_id(long size, double alpha, double dens, uint32_t seed);

// Read the settings from a job name; false (with the reason in `error`) if it is not a
// job name written as the sweep writes it.
bool parse_job_id(const std::string& id, Job& job, std::string& error);

// The "engine" field of a written file: which native build wrote it, and in which mode
// (`name`, "lgca_native run" for run.json).
std::string engine_json(const char* name);

// Append `item` to a comma-separated list.
void append_item(std::string& list, const std::string& item);

// Whether s is `length` lowercase hexadecimal digits.
bool is_lower_hex(const std::string& s, std::size_t length);

// Whether v is a number with an integer value from min_value to max_value.
bool is_integer(const JsonValue* v, double min_value, double max_value);

}  // namespace lgca
