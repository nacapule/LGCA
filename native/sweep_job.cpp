// native/sweep_job.cpp — what the `run` and `trace` modes share. See sweep_job.hpp.

#include "sweep_job.hpp"

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "measure.hpp"

// Set by build.sh: the SHA-256 of the native sources and the compiler flags, recorded in
// the "engine" field. A build without them records "unknown".
#ifndef LGCA_SOURCES_SHA256
#define LGCA_SOURCES_SHA256 "unknown"
#endif
#ifndef LGCA_BUILD_FLAGS
#define LGCA_BUILD_FLAGS "unknown"
#endif

namespace lgca {

std::vector<long> sample_steps(long t) {
    std::vector<long> steps;
    for (long s = 0; s <= t; s += s < 5000 ? 25 : 100) steps.push_back(s);
    return steps;
}

std::vector<long> tick_steps(long t) {
    std::vector<long> steps;
    for (const long tick : TICKS)
        if (tick <= t) steps.push_back(tick);
    return steps;
}

std::string state_file_name(long t) { return "state-t" + std::to_string(t) + ".bin.gz"; }

bool is_state_file_name(const std::string& name, long* t) {
    const std::string head = "state-t", tail = ".bin.gz";
    if (name.size() <= head.size() + tail.size() || name.compare(0, head.size(), head) != 0 ||
        name.compare(name.size() - tail.size(), tail.size(), tail) != 0)
        return false;
    for (std::size_t p = head.size(); p < name.size() - tail.size(); ++p)
        if (name[p] < '0' || name[p] > '9') return false;
    if (t) {
        const std::string digits = name.substr(head.size(), name.size() - head.size() - tail.size());
        long long value = 0;
        *t = digits.size() <= 18 && parse_integer(digits, 0, 1000000000000000000LL, value)
                 ? static_cast<long>(value)
                 : -1;
    }
    return true;
}

Params job_params(const Job& job) {
    Params params;
    params.W = params.H = static_cast<int>(job.size);
    params.model = Model::Boson;
    params.dens = job.dens;
    params.seed = job.seed;
    params.kernel = Kernel::Power;
    params.alpha = job.alpha;
    params.boson_field = BosonField::Both;
    params.boson_align = BosonAlign::Polar;
    return params;
}

// The whole text must be a decimal integer from min_value to max_value (strtoll alone
// would also accept leading spaces and a '+').
bool parse_integer(const std::string& text, long long min_value, long long max_value,
                   long long& value) {
    if (text.empty() || !(std::isdigit(static_cast<unsigned char>(text[0])) || text[0] == '-'))
        return false;
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(text.c_str(), &end, 10);
    if (errno != 0 || *end != '\0' || v < min_value || v > max_value) return false;
    value = v;
    return true;
}

// The whole text must be a finite number (as for --alpha in lgca_cli.cpp).
bool parse_number(const std::string& text, double& value) {
    if (text.empty() || std::isspace(static_cast<unsigned char>(text[0]))) return false;
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(text.c_str(), &end);
    if ((errno != 0 && errno != ERANGE) || *end != '\0' || !std::isfinite(v)) return false;
    value = v;
    return true;
}

// The job name the sweep uses: both-L<size>-a<alpha>-d<dens>-seed<seed>, each number
// written as JavaScript writes it in a template string (the same text as json_number).
std::string job_id(long size, double alpha, double dens, uint32_t seed) {
    return "both-L" + std::to_string(size) + "-a" + json_number(alpha) + "-d" + json_number(dens) +
           "-seed" + std::to_string(seed);
}

// Read the settings from a job name. The JavaScript runner matches
// /^both-L(\d+)-a([\d.]+)-d([\d.]+)-seed(\d+)$/ and converts each part with Number().
// Here a part must also be written as JavaScript writes that number (no leading zeros,
// no trailing ".0"), so that the name maps to one job and back to the same name.
bool parse_job_id(const std::string& id, Job& job, std::string& error) {
    error = "job name must look like both-L90-a0.85-d0.4-seed12345";
    const auto take = [&id](std::size_t& p, const char* word) {
        const std::size_t n = std::strlen(word);
        if (id.compare(p, n, word) != 0) return false;
        p += n;
        return true;
    };
    // The longest run of characters from `allowed` starting at p (at least one).
    const auto part = [&id](std::size_t& p, const char* allowed, std::string& out) {
        const std::size_t from = p;
        while (p < id.size() && id[p] != '\0' && std::strchr(allowed, id[p])) ++p;
        out = id.substr(from, p - from);
        return p > from;
    };
    std::size_t p = 0;
    std::string size_text, alpha_text, dens_text, seed_text;
    if (!take(p, "both-L") || !part(p, "0123456789", size_text) || !take(p, "-a") ||
        !part(p, "0123456789.", alpha_text) || !take(p, "-d") ||
        !part(p, "0123456789.", dens_text) || !take(p, "-seed") ||
        !part(p, "0123456789", seed_text) || p != id.size())
        return false;
    long long size, seed;
    double alpha, dens;
    if (!parse_integer(size_text, 0, MAX_SITES, size) ||
        !parse_integer(seed_text, 0, 4294967295LL, seed) || !parse_number(alpha_text, alpha) ||
        !parse_number(dens_text, dens))
        return false;
    job.id = id;
    job.size = static_cast<long>(size);
    job.alpha = alpha;
    job.dens = dens;
    job.seed = static_cast<uint32_t>(seed);
    if (job_id(job.size, job.alpha, job.dens, job.seed) != id) {
        error = "job name " + id + " is not written as the sweep writes it (it would be " +
                job_id(job.size, job.alpha, job.dens, job.seed) + ")";
        return false;
    }
    return true;
}

// The "engine" field: which native build wrote the file.
std::string engine_json(const char* name) {
    return std::string("{\"name\":") + json_string(name) + ",\"sourcesSha256\":" +
           json_string(LGCA_SOURCES_SHA256) + ",\"compiler\":" + json_string(__VERSION__) +
           ",\"flags\":" + json_string(LGCA_BUILD_FLAGS) + "}";
}

// Append `item` to a comma-separated list.
void append_item(std::string& list, const std::string& item) {
    if (!list.empty()) list += ',';
    list += item;
}

bool is_lower_hex(const std::string& s, std::size_t length) {
    if (s.size() != length) return false;
    for (const char c : s)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

// Whether v is a number with an integer value from min_value to max_value.
bool is_integer(const JsonValue* v, double min_value, double max_value) {
    return v && v->is_number() && v->number >= min_value && v->number <= max_value &&
           v->number == std::floor(v->number);
}

}  // namespace lgca
