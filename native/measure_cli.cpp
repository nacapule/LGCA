// native/measure_cli.cpp — the native program's `measure` mode.
//
// Usage:
//
//   lgca_native measure [--W N --H N] [--factor F] [--snap | --island] [--t T] [--repeat R] FILE...
//
// Measures saved lattices with the sweeps' observables (measure.hpp: the same numbers
// as the JavaScript measurements() function, bit for bit) and
// prints one line of JSON per FILE, in the shapes of the sweep's run.json:
//
//   default   a `series` row   {"t":T,"N":...,"mean":...,...,"clusterPolar":...}
//   --snap    a `snaps` entry  {"t":T,"file":...,"sha256":...,"rng":[...],
//                               "N":...,...,"clusterPolar":...,"thresholdChecks":[...]}
//   --island  the island (measure.hpp, Island)  {"t":T,"file":...,"N":...,"nmax":...,
//                               "seed":...,"islandMass":...,...,"largestMass":...,
//                               "largestArea":...,...,"densest2Area":...}
//
// The text is what JSON.stringify writes for the same values, character for character.
//
// FILE is one of:
//   - a raw lattice: the W*H*6 int32 values occ[(j*W + i)*6 + k], little-endian, as the
//     sweeps save them (state-t<T>.bin.gz in each job's folder);
//   - a state file of this program (--dump-state): 33 generator words, then the lattice.
//     Its generator words become the snap's "rng" field; a raw lattice has none, and
//     its snap has no "rng" field.
// Either may be gzip-compressed: a file that starts with the gzip bytes 1f 8b 08 is
// decompressed as Node's gunzipSync does it (several gzip members in a row are read as
// one; after a member, a zero byte ends the reading and the rest of the file is ignored,
// whatever it holds; any other byte must begin another member; a damaged or cut-off
// stream is an error). Sizes are capped: the content at the size --W and --H give (else
// 10^8 sites), and the gzip bytes read at slightly more than that, so that endless input
// stops with an error (file_io.hpp, read_file). A plain state file whose first generator
// word happens to start with those 3 bytes (1 file in 16.7 million) is refused as a
// damaged gzip stream; compress it with gzip to measure it.
//
// Options:
//   --W N, --H N   the lattice size. Without them the lattice is taken to be square and
//                  its side is found from the file size.
//   --factor F     the cluster threshold is max(1, ceil(F * N / (W*H))); default 2, the
//                  sweeps' value. The snaps' thresholdChecks always use 1.5 and 4.
//   --snap         write snaps entries instead of series rows. "file" is FILE's name
//                  without its folder, "sha256" the SHA-256 of the raw lattice bytes.
//   --island       write the island instead: the busiest site and every occupied site
//                  connected to it (measure.hpp). Not with --snap. A lattice with a
//                  negative channel count is refused.
//   --t T          the step count written as "t". Without it, t is read from the file
//                  name: the number after the last ".t" or "-t" (state-t500.bin.gz,
//                  run.r0.t500).
//   --repeat R     measure each lattice R times and print, on stderr, the median time
//                  of one measurement (and of one whole snap with --snap, of the
//                  island alone with --island): the measurement cost per sample. The
//                  JSON on stdout is the same.
// Options may also be written --option=value. Errors stop the program with a message
// and exit status 1.

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "file_io.hpp"
#include "measure.hpp"
#include "sha256.hpp"

namespace lgca {

namespace {

constexpr std::size_t RNG_BYTES = 33 * sizeof(uint32_t);   // a state file's generator words
constexpr std::size_t SITE_BYTES = 6 * sizeof(int32_t);    // one site's 6 channels

bool parse_long(const std::string& text, long min_value, long& value) {
    if (text.empty()) return false;
    char* end = nullptr;
    errno = 0;
    const long v = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || *end != '\0' || v < min_value || v > 1000000000L) return false;
    value = v;
    return true;
}

// Whether a file name is valid UTF-8 (no overlong forms, surrogates or values past
// U+10FFFF). A snap writes the name as a JSON string, which must be text.
bool valid_utf8(const std::string& s) {
    std::size_t p = 0;
    while (p < s.size()) {
        const unsigned char c = static_cast<unsigned char>(s[p]);
        int extra;
        uint32_t value;
        if (c < 0x80) {
            ++p;
            continue;
        } else if (c >= 0xc2 && c <= 0xdf) {
            extra = 1;
            value = c & 0x1f;
        } else if (c >= 0xe0 && c <= 0xef) {
            extra = 2;
            value = c & 0x0f;
        } else if (c >= 0xf0 && c <= 0xf4) {
            extra = 3;
            value = c & 0x07;
        } else {
            return false;
        }
        if (s.size() - p - 1 < static_cast<std::size_t>(extra)) return false;
        for (int q = 1; q <= extra; ++q) {
            const unsigned char d = static_cast<unsigned char>(s[p + q]);
            if ((d & 0xc0) != 0x80) return false;
            value = (value << 6) | (d & 0x3f);
        }
        if ((extra == 2 && value < 0x800) || (extra == 3 && value < 0x10000) ||
            (value >= 0xd800 && value <= 0xdfff) || value > 0x10ffff) {
            return false;
        }
        p += 1 + extra;
    }
    return true;
}

// t from a file name: the digits after the last ".t" or "-t" that run to the end of
// the name or to a ".".
std::optional<long> step_from_name(const std::string& name) {
    std::optional<long> found;
    for (std::size_t p = 1; p < name.size(); ++p) {
        if (name[p] != 't' || (name[p - 1] != '.' && name[p - 1] != '-')) continue;
        std::size_t q = p + 1;
        while (q < name.size() && name[q] >= '0' && name[q] <= '9') ++q;
        long t = 0;
        if (q > p + 1 && (q == name.size() || name[q] == '.') &&
            parse_long(name.substr(p + 1, q - p - 1), 0, t)) {
            found = t;
        }
    }
    return found;
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size() / 2;
    return v.size() % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

void usage() {
    std::cerr << "usage: lgca_native measure [--W N --H N] [--factor F] [--snap | --island]"
                 " [--t T] [--repeat R] FILE...\n";
}

}  // namespace

int measure_main(int argc, char** argv) {
    long W = 0, H = 0, fixed_t = -1, repeat = 1;
    double factor = 2;
    bool snap = false, island = false;
    std::vector<std::string> files;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) {
            files.push_back(arg);
            continue;
        }
        std::string name = arg.substr(2), value;
        if (name == "snap") {
            snap = true;
            continue;
        }
        if (name == "island") {
            island = true;
            continue;
        }
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
        if (name == "snap" || name == "island") {
            std::cerr << "error: --" << name << " takes no value\n";
            return 1;
        } else if (name == "W" || name == "H") {
            if (!parse_long(value, 1, name == "W" ? W : H)) {
                std::cerr << "error: --" << name << " must be a positive integer\n";
                return 1;
            }
        } else if (name == "factor") {
            char* end = nullptr;
            factor = std::strtod(value.c_str(), &end);
            if (value.empty() || *end != '\0' || !std::isfinite(factor)) {
                std::cerr << "error: --factor must be a finite number\n";
                return 1;
            }
        } else if (name == "t") {
            if (!parse_long(value, 0, fixed_t)) {
                std::cerr << "error: --t must be an integer from 0 to 1000000000\n";
                return 1;
            }
        } else if (name == "repeat") {
            if (!parse_long(value, 1, repeat)) {
                std::cerr << "error: --repeat must be a positive integer\n";
                return 1;
            }
        } else {
            std::cerr << "error: unknown option --" << name << "\n";
            usage();
            return 1;
        }
    }
    if (files.empty()) {
        usage();
        return 1;
    }
    if (snap && island) {
        std::cerr << "error: give --snap or --island, not both\n";
        return 1;
    }
    if ((W == 0) != (H == 0)) {
        std::cerr << "error: give both --W and --H, or neither (a square lattice)\n";
        return 1;
    }
    if (W > 0 && W * H > 100000000L) {
        std::cerr << "error: W*H must be at most 100000000\n";
        return 1;
    }
    // The saved lattices are little-endian; so is every machine this program targets.
    const uint32_t one = 1;
    if (*reinterpret_cast<const unsigned char*>(&one) != 1) {
        std::cerr << "error: this program reads little-endian lattices and needs a"
                     " little-endian machine\n";
        return 1;
    }

    std::optional<Measurer> measurer;
    std::vector<unsigned char> bytes;
    // The most bytes a file may hold: a state file of the given size, or of the largest
    // square lattice.
    const std::size_t limit =
        static_cast<std::size_t>(W > 0 ? W * H : 100000000L) * SITE_BYTES + RNG_BYTES;
    for (const std::string& path : files) {
        std::string error;
        if (!read_file(path, limit, bytes, error)) {
            std::cerr << "error: " << path << ": cannot read: " << error << "\n";
            return 1;
        }
        // A raw lattice is a multiple of 24 bytes; a state file is 132 bytes more, and
        // 132 is not a multiple of 24, so the size tells the two apart.
        const bool has_rng = bytes.size() >= RNG_BYTES &&
                             bytes.size() % SITE_BYTES == RNG_BYTES % SITE_BYTES;
        const std::size_t lattice_bytes = bytes.size() - (has_rng ? RNG_BYTES : 0);
        long w = W, h = H;
        if (w == 0) {   // square: side = sqrt(sites)
            const std::size_t sites = lattice_bytes / SITE_BYTES;
            const long side = std::lround(std::sqrt(static_cast<double>(sites)));
            if (lattice_bytes % SITE_BYTES != 0 || sites == 0 || sites > 100000000 ||
                static_cast<std::size_t>(side) * side != sites) {
                std::cerr << "error: " << path << ": " << bytes.size() << " bytes is not a square"
                             " lattice (or a state file of one); give --W and --H\n";
                return 1;
            }
            w = h = side;
        }
        const std::size_t expected = static_cast<std::size_t>(w) * h * SITE_BYTES;
        if (lattice_bytes != expected) {
            std::cerr << "error: " << path << ": " << bytes.size() << " bytes; a " << w << " x " << h
                      << " lattice has " << expected << " bytes, a state file "
                      << expected + RNG_BYTES << "\n";
            return 1;
        }
        const std::string name = path.substr(path.find_last_of('/') + 1);
        if ((snap || island) && !valid_utf8(name)) {
            std::cerr << "error: " << path << ": the file name is not UTF-8 text, and the"
                         " output writes it as a JSON string\n";
            return 1;
        }
        long t = fixed_t;
        if (t < 0) {
            const std::optional<long> from_name = step_from_name(name);
            if (!from_name) {
                std::cerr << "error: " << path << ": no step count in the file name"
                             " (like state-t500.bin.gz); give --t T\n";
                return 1;
            }
            t = *from_name;
        }

        // The lattice as int32 values (copied: the byte buffer need not be aligned).
        std::vector<int32_t> occ(lattice_bytes / sizeof(int32_t));
        std::memcpy(occ.data(), bytes.data() + (has_rng ? RNG_BYTES : 0), lattice_bytes);
        uint32_t rng[33];
        if (has_rng) std::memcpy(rng, bytes.data(), RNG_BYTES);

        if (!measurer || measurer->width() != w || measurer->height() != h) {
            measurer.emplace(static_cast<int>(w), static_cast<int>(h));
        }
        using Clock = std::chrono::steady_clock;
        std::vector<double> series_seconds, snap_seconds, island_seconds;
        std::string line;
        for (long r = 0; r < repeat; ++r) {
            const Clock::time_point start = Clock::now();
            const Measurements m = measurer->measure(occ.data(), factor);
            const Clock::time_point measured = Clock::now();
            if (island) {
                Island s;
                try {
                    s = measurer->island(occ.data());
                } catch (const std::exception& e) {
                    std::cerr << "error: " << path << ": " << e.what() << "\n";
                    return 1;
                }
                island_seconds.push_back(
                    std::chrono::duration<double>(Clock::now() - measured).count());
                line = "{\"t\":" + json_number(static_cast<double>(t)) + ",\"file\":" +
                       json_string(name) + ",\"N\":" + json_number(m.N) +
                       ",\"nmax\":" + json_number(m.nmax) + island_fields_json(s) + "}";
            } else if (snap) {
                const std::vector<ThresholdCheck> checks = measurer->threshold_checks(occ.data());
                const std::string hash = sha256_hex(occ.data(), lattice_bytes);
                line = snap_json(t, name, hash, has_rng ? rng : nullptr, m, checks);
            } else {
                line = series_json(t, m);
            }
            const Clock::time_point done = Clock::now();
            series_seconds.push_back(std::chrono::duration<double>(measured - start).count());
            snap_seconds.push_back(std::chrono::duration<double>(done - start).count());
        }
        std::fwrite(line.data(), 1, line.size(), stdout);
        std::fputc('\n', stdout);
        if (repeat > 1) {
            std::fprintf(stderr, "%s: %ld x %ld sites, median of %ld: %.6g s per measurement%s",
                         name.c_str(), w, h, repeat, median(series_seconds),
                         snap || island ? "" : "\n");
            if (snap) std::fprintf(stderr, ", %.6g s per snap\n", median(snap_seconds));
            if (island) std::fprintf(stderr, ", %.6g s per island\n", median(island_seconds));
        }
    }
    if (std::fflush(stdout) != 0) {
        std::cerr << "error: cannot write the output\n";
        return 1;
    }
    return 0;
}

}  // namespace lgca
