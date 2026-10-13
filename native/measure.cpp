// native/measure.cpp — the sweeps' per-sample observables. See measure.hpp.
//
// Each step below quotes the line of the JavaScript measurements() function (see
// measure.hpp) that it reproduces. Keep the order of every sum as it is: floating-point
// addition is not associative, and a reordered sum would change the last bits of the
// results.

#include "measure.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace lgca {

namespace {

// The 6 neighbour offsets (di, dj) in the order the cluster search visits them
// (NB in the JavaScript function).
constexpr int NB[6][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, -1}, {-1, 1}};

// n[x] += v on an Int32Array: the sum taken modulo 2^32 and read back as a signed 32-bit
// integer. (Unsigned arithmetic in C++ wraps the same way; signed overflow would be
// undefined behaviour.)
inline int32_t int32_add(int32_t a, int32_t b) {
    return static_cast<int32_t>(static_cast<uint32_t>(a) + static_cast<uint32_t>(b));
}

// The series fields, "N":... to "clusterPolar":..., each preceded by a comma.
std::string fields_json(const Measurements& m) {
    std::string s;
    const auto field = [&s](const char* name, double v) {
        s += ",\"";
        s += name;
        s += "\":";
        s += json_number(v);
    };
    field("N", m.N);
    field("mean", m.mean);
    field("nmax", m.nmax);
    field("kmax", m.kmax);
    field("siteFraction", m.siteFraction);
    field("channelFraction", m.channelFraction);
    field("effectiveArea", m.effectiveArea);
    field("occupied", static_cast<double>(m.occupied));
    field("polar", m.polar);
    field("localPolar", m.localPolar);
    field("weightedLocalPolar", m.weightedLocalPolar);
    field("threshold", m.threshold);
    field("clusters", static_cast<double>(m.clusters));
    field("clusterMass", m.clusterMass);
    field("clusterArea", static_cast<double>(m.clusterArea));
    field("clusterDensity", m.clusterDensity);
    field("clusterFraction", m.clusterFraction);
    field("clusterPolar", m.clusterPolar);
    return s;
}

}  // namespace

std::string json_string(const std::string& s) {
    std::string out = "\"";
    for (const unsigned char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char code[8];
                    std::snprintf(code, sizeof code, "\\u%04x", c);
                    out += code;
                } else {
                    out += static_cast<char>(c);
                }
        }
    }
    return out + "\"";
}

Measurer::Measurer(int W, int H) : W_(W), H_(H) {
    // Site indices are kept in 32-bit integers, as in the JavaScript function (an
    // Int32Array queue).
    if (W < 1 || H < 1 || static_cast<int64_t>(W) * H > std::numeric_limits<int32_t>::max()) {
        throw std::invalid_argument("Measurer: W and H must be positive and W*H at most 2^31 - 1");
    }
    const std::size_t A = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
    n_.resize(A);
    px_.resize(A);
    py_.resize(A);
    seen_.resize(A);
    queue_.resize(A);
    best_.reserve(A);
}

Measurements Measurer::measure(const int32_t* occ, double factor) {
    const int32_t A = W_ * H_;
    // Channel velocities: cx=[1,.5,-.5,-1,-.5,.5], cy=[0,Math.sqrt(3)/2,...]
    const double s = std::sqrt(3.0) / 2;
    const double cx[6] = {1, .5, -.5, -1, -.5, .5};
    const double cy[6] = {0, s, s, 0, -s, -s};

    // ---- per site: n, px, py, and the lattice totals ---------------------------------
    // let N=0,nmax=0,kmax=0,sum2=0,occupied=0,lpol=0,weighted=0,gx=0,gy=0
    double N = 0, sum2 = 0, lpol = 0, weighted = 0, gx = 0, gy = 0;
    int32_t nmax = 0, kmax = 0;
    int64_t occupied = 0;
    for (int32_t x = 0; x < A; ++x) {
        int32_t n = 0;
        double px = 0, py = 0;
        // for(let k=0;k<6;k++){const v=occ[6*x+k];n[x]+=v;px[x]+=v*cx[k];py[x]+=v*cy[k];kmax=Math.max(kmax,v);}
        for (int k = 0; k < 6; ++k) {
            const int32_t v = occ[6 * static_cast<int64_t>(x) + k];
            n = int32_add(n, v);
            px += static_cast<double>(v) * cx[k];
            py += static_cast<double>(v) * cy[k];
            kmax = std::max(kmax, v);
        }
        n_[x] = n;
        px_[x] = px;
        py_[x] = py;
        // N+=n[x];sum2+=n[x]*n[x];nmax=Math.max(nmax,n[x]);
        N += n;
        sum2 += static_cast<double>(n) * static_cast<double>(n);   // a double product, as in JavaScript
        nmax = std::max(nmax, n);
        // const p=Math.hypot(px[x],py[x]);gx+=px[x];gy+=py[x];weighted+=p;if(n[x]){occupied++;lpol+=p/n[x];}
        const double p = js_hypot(px, py);
        gx += px;
        gy += py;
        weighted += p;
        if (n != 0) {
            ++occupied;
            lpol += p / n;
        }
    }

    // ---- the threshold: Math.max(1, Math.ceil(thresholdMultiplier*N/A)) --------------
    // Math.max returns NaN if an argument is NaN; std::max would return 1. (The ceiling
    // is NaN only for an infinite or NaN factor, which the command line refuses.)
    const double ceiling = std::ceil(factor * N / A);
    const double threshold = std::isnan(ceiling) ? ceiling : std::max(1.0, ceiling);

    // ---- clusters: breadth-first search from every unvisited site with n >= threshold --
    // neigh(x,d) = ((floor(x/W)+NB[d][1]+H)%H)*W + (x%W+NB[d][0]+W)%W: the neighbour of
    // site (i, j) in direction d, with periodic boundaries.
    std::fill(seen_.begin(), seen_.end(), 0);
    double mass = 0;
    int64_t clusters = 0;
    best_.clear();
    for (int32_t x = 0; x < A; ++x) {
        if (seen_[x] || !(n_[x] >= threshold)) continue;
        int32_t head = 0, tail = 1;
        double m = 0;
        queue_[0] = x;
        seen_[x] = 1;
        while (head < tail) {
            const int32_t v = queue_[head++];
            m += n_[v];
            const int32_t i = v % W_, j = v / W_;
            for (int d = 0; d < 6; ++d) {
                // The same site as neigh(v,d): i + di and j + dj wrapped around the
                // lattice (they are at most one step outside it).
                int32_t zi = i + NB[d][0], zj = j + NB[d][1];
                if (zi < 0) zi += W_; else if (zi >= W_) zi -= W_;
                if (zj < 0) zj += H_; else if (zj >= H_) zj -= H_;
                const int32_t z = zj * W_ + zi;
                if (!seen_[z] && n_[z] >= threshold) {
                    seen_[z] = 1;
                    queue_[tail++] = z;
                }
            }
        }
        ++clusters;
        // if(m>mass){mass=m;best=Array.from(queue.subarray(0,t));}  (ties: the first wins)
        if (m > mass) {
            mass = m;
            best_.assign(queue_.begin(), queue_.begin() + tail);
        }
    }
    // let bx=0,by=0;for(const x of best){bx+=px[x];by+=py[x];}
    double bx = 0, by = 0;
    for (const int32_t x : best_) {
        bx += px_[x];
        by += py_[x];
    }

    // ---- the result object, field by field ------------------------------------------
    // (a JavaScript test `N ? ... : 0` is false for 0 only here: these totals are never NaN)
    const double area = static_cast<double>(best_.size());
    Measurements r;
    r.N = N;
    r.mean = N / A;
    r.nmax = nmax;
    r.kmax = kmax;
    r.siteFraction = N != 0 ? nmax / N : 0;
    r.channelFraction = N != 0 ? kmax / N : 0;
    r.effectiveArea = sum2 != 0 ? N * N / sum2 : 0;
    r.occupied = occupied;
    r.polar = N != 0 ? js_hypot(gx, gy) / N : 0;
    r.localPolar = occupied != 0 ? lpol / static_cast<double>(occupied) : 0;
    r.weightedLocalPolar = N != 0 ? weighted / N : 0;
    r.threshold = threshold;
    r.clusters = clusters;
    r.clusterMass = mass;
    r.clusterArea = static_cast<int64_t>(best_.size());
    r.clusterDensity = area != 0 ? mass / area : 0;
    r.clusterFraction = N != 0 ? mass / N : 0;
    r.clusterPolar = mass != 0 ? js_hypot(bx, by) / mass : 0;
    return r;
}

std::vector<ThresholdCheck> Measurer::threshold_checks(const int32_t* occ) {
    // [1.5,4].map(f=>{const c=measurements(sim.occ,size,size,f);
    //   return {factor:f,threshold:c.threshold,mass:c.clusterMass,area:c.clusterArea,density:c.clusterDensity}})
    std::vector<ThresholdCheck> checks;
    for (const double f : THRESHOLD_CHECK_FACTORS) {
        const Measurements c = measure(occ, f);
        checks.push_back({f, c.threshold, c.clusterMass, c.clusterArea, c.clusterDensity});
    }
    return checks;
}

Island Measurer::island(const int32_t* occ) {
    const int32_t A = W_ * H_;
    const double s = std::sqrt(3.0) / 2;
    const double cx[6] = {1, .5, -.5, -1, -.5, .5};
    const double cy[6] = {0, s, s, 0, -s, -s};
    if (ui_.empty()) {
        ui_.resize(static_cast<std::size_t>(A));
        uj_.resize(static_cast<std::size_t>(A));
    }

    // Particles and momentum per site; the busiest site (the first of equal ones).
    // Unlike measure(), which reproduces JavaScript's wrapping Int32Array sums, a site
    // total that does not fit in 32 bits is refused rather than wrapped.
    Island r;
    int32_t nmax = 0;
    for (int32_t x = 0; x < A; ++x) {
        int64_t n = 0;
        double px = 0, py = 0;
        for (int k = 0; k < 6; ++k) {
            const int32_t v = occ[6 * static_cast<int64_t>(x) + k];
            // A negative channel could hide inside a positive site total; no lattice of
            // the engine has one, so it is refused (measure() keeps JavaScript's sums).
            if (v < 0) throw std::invalid_argument("island: a channel count is negative");
            n += v;
            px += static_cast<double>(v) * cx[k];
            py += static_cast<double>(v) * cy[k];
        }
        if (n > std::numeric_limits<int32_t>::max())
            throw std::overflow_error("island: a site total is above 2^31 - 1");
        n_[x] = static_cast<int32_t>(n);
        px_[x] = px;
        py_[x] = py;
        if (n > nmax) {
            nmax = n;
            r.seed = x;
        }
    }
    if (r.seed < 0) return r;   // no particles, no island

    // The neighbour of (i, j) in direction d, wrapped into the box.
    const auto neighbour = [this](int32_t i, int32_t j, int d) {
        int32_t zi = i + NB[d][0], zj = j + NB[d][1];
        if (zi < 0) zi += W_; else if (zi >= W_) zi -= W_;
        if (zj < 0) zj += H_; else if (zj >= H_) zj -= H_;
        return zj * W_ + zi;
    };

    // The island: breadth-first search over occupied sites from the busiest one,
    // carrying unwrapped coordinates (ui, uj).
    std::fill(seen_.begin(), seen_.end(), 0);
    int32_t head = 0, tail = 1;
    queue_[0] = r.seed;
    seen_[r.seed] = 1;
    ui_[r.seed] = r.seed % W_;
    uj_[r.seed] = r.seed / W_;
    double bx = 0, by = 0;
    int32_t island_low = r.seed;   // the island's lowest site index
    while (head < tail) {
        const int32_t v = queue_[head++];
        r.mass += n_[v];
        bx += px_[v];
        by += py_[v];
        island_low = std::min(island_low, v);
        const int32_t i = v % W_, j = v / W_;
        for (int d = 0; d < 6; ++d) {
            const int32_t z = neighbour(i, j, d);
            if (n_[z] == 0) continue;
            const int32_t zi = ui_[v] + NB[d][0], zj = uj_[v] + NB[d][1];
            if (!seen_[z]) {
                seen_[z] = 1;
                ui_[z] = zi;
                uj_[z] = zj;
                queue_[tail++] = z;
            } else if (ui_[z] != zi || uj_[z] != zj) {
                r.wraps = true;   // reached again by a path once around the box
            }
        }
    }
    r.area = tail;
    r.polar = r.mass != 0 ? js_hypot(bx, by) / r.mass : 0;
    if (r.wraps) {
        r.ci = r.cj = r.rg = std::numeric_limits<double>::quiet_NaN();
    } else {
        // Centre and radius of gyration in Cartesian coordinates: site (i, j) sits at
        // i*(1, 0) + j*(1/2, -sqrt(3)/2), the streaming offsets' geometry (PHYSICS.md).
        double si = 0, sj = 0;
        for (int32_t q = 0; q < tail; ++q) {
            const int32_t v = queue_[q];
            si += static_cast<double>(n_[v]) * ui_[v];
            sj += static_cast<double>(n_[v]) * uj_[v];
        }
        const double mi = si / r.mass, mj = sj / r.mass;
        double r2 = 0;
        for (int32_t q = 0; q < tail; ++q) {
            const int32_t v = queue_[q];
            const double di = ui_[v] - mi, dj = uj_[v] - mj;
            const double x = di + 0.5 * dj, y = -s * dj;
            r2 += static_cast<double>(n_[v]) * (x * x + y * y);
        }
        r.rg = std::sqrt(r2 / r.mass);
        r.ci = std::fmod(std::fmod(mi, W_) + W_, W_);
        r.cj = std::fmod(std::fmod(mj, H_) + H_, H_);
    }

    // The heaviest, the densest and the densest of at least two sites among all islands.
    // Islands are taken in the order of their lowest site index, and a choice changes only
    // for a larger value or for an equal value of an island earlier in that order: the
    // first of equal values in heaviest.py's order (measure.hpp).
    struct Choice {
        bool any = false;
        double value = 0;   // mass, or mass / area
        int32_t low = 0;    // the island's lowest site index
        double mass = 0;
        int64_t area = 0;
    };
    Choice heaviest, densest, densest2;
    const auto consider = [](Choice& c, double value, int32_t low, double mass, int64_t area) {
        if (!c.any || value > c.value || (value == c.value && low < c.low))
            c = {true, value, low, mass, area};
    };
    const auto rank = [&](int32_t low, double mass, int64_t area) {
        const double density = mass / static_cast<double>(area);
        consider(heaviest, mass, low, mass, area);
        consider(densest, density, low, mass, area);
        if (area >= 2) consider(densest2, density, low, mass, area);
    };
    rank(island_low, r.mass, r.area);

    // Every island: the same search from each occupied site not yet reached.
    r.components = 1;
    for (int32_t x = 0; x < A; ++x) {
        if (seen_[x] || n_[x] == 0) continue;
        head = 0;
        tail = 1;
        queue_[0] = x;
        seen_[x] = 1;
        double m = 0;
        while (head < tail) {
            const int32_t v = queue_[head++];
            m += n_[v];
            const int32_t i = v % W_, j = v / W_;
            for (int d = 0; d < 6; ++d) {
                const int32_t z = neighbour(i, j, d);
                if (!seen_[z] && n_[z] != 0) {
                    seen_[z] = 1;
                    queue_[tail++] = z;
                }
            }
        }
        ++r.components;
        // The scan meets each of these islands first at its lowest site, x. (The busiest
        // site's island, searched before the scan, kept its lowest site in island_low.)
        rank(x, m, tail);
    }
    r.largestMass = heaviest.mass;
    r.largestArea = heaviest.area;
    r.densestMass = densest.mass;
    r.densestArea = densest.area;
    if (densest2.any) {
        r.densest2Mass = densest2.mass;
        r.densest2Area = densest2.area;
    }
    return r;
}

std::string island_fields_json(const Island& s) {
    std::string out;
    const auto field = [&out](const char* name, double v) {
        out += ",\"";
        out += name;
        out += "\":";
        out += json_number(v);
    };
    field("seed", s.seed);
    field("islandMass", s.mass);
    field("islandArea", static_cast<double>(s.area));
    field("islandPolar", s.polar);
    field("wraps", s.wraps ? 1 : 0);
    field("ci", s.ci);
    field("cj", s.cj);
    field("rg", s.rg);
    field("components", static_cast<double>(s.components));
    field("largestMass", s.largestMass);
    field("largestArea", static_cast<double>(s.largestArea));
    field("densestMass", s.densestMass);
    field("densestArea", static_cast<double>(s.densestArea));
    field("densest2Mass", s.densest2Mass);
    field("densest2Area", island_densest2_area(s));
    return out;
}

double js_hypot(double a, double b) {
    // V8's MathHypot (src/builtins/math.tq) for two arguments, step by step.
    const double values[2] = {a, b};
    double abs_values[2] = {0, 0};
    bool one_arg_is_nan = false;
    double max = 0;
    for (int i = 0; i < 2; ++i) {
        if (std::isnan(values[i])) {
            one_arg_is_nan = true;
        } else {
            abs_values[i] = std::fabs(values[i]);
            if (abs_values[i] > max) max = abs_values[i];
        }
    }
    if (max == std::numeric_limits<double>::infinity()) return max;
    if (one_arg_is_nan) return std::numeric_limits<double>::quiet_NaN();
    if (max == 0) return 0;
    // Kahan summation of the squares, each value normalized by the largest one.
    double sum = 0, compensation = 0;
    for (int i = 0; i < 2; ++i) {
        const double n = abs_values[i] / max;
        const double summand = n * n - compensation;
        const double preliminary = sum + summand;
        compensation = (preliminary - sum) - summand;
        sum = preliminary;
    }
    return std::sqrt(sum) * max;
}

std::string json_number(double v) {
    if (!std::isfinite(v)) return "null";   // JSON.stringify(NaN) and of ±Infinity
    if (v == 0) return "0";                 // 0 and -0
    // The shortest digits that read back as v; among equally short ones the closest to
    // v (std::to_chars guarantees both, as ECMAScript's Number::toString requires).
    char buf[64];
    const std::to_chars_result res =
        std::to_chars(buf, buf + sizeof buf, v, std::chars_format::scientific);
    const std::string sci(buf, res.ptr);   // e.g. "-1.2345e-07"
    const bool negative = sci[0] == '-';
    const std::size_t first = negative ? 1 : 0, e = sci.find('e');
    std::string digits = sci.substr(first, 1);
    if (e > first + 1) digits += sci.substr(first + 2, e - first - 2);
    const int k = static_cast<int>(digits.size());
    // v = 0.d1d2...dk * 10^n
    const int n = std::atoi(sci.c_str() + e + 1) + 1;

    // ECMAScript Number::toString(v), for v > 0 (the sign is added back at the end).
    std::string out;
    if (k <= n && n <= 21) {               // an integer: the digits, then n - k zeros
        out = digits + std::string(n - k, '0');
    } else if (0 < n && n <= 21) {         // a decimal point inside the digits
        out = digits.substr(0, n) + "." + digits.substr(n);
    } else if (-6 < n && n <= 0) {         // 0.000ddd
        out = "0." + std::string(-n, '0') + digits;
    } else {                               // exponent form: d.ddde+X or de-X
        out = digits.substr(0, 1);
        if (k > 1) out += "." + digits.substr(1);
        out += n - 1 >= 0 ? "e+" : "e-";
        out += std::to_string(std::abs(n - 1));
    }
    return negative ? "-" + out : out;
}

std::string series_json(long t, const Measurements& m) {
    return "{\"t\":" + json_number(static_cast<double>(t)) + fields_json(m) + "}";
}

std::string snap_json(long t, const std::string& file, const std::string& sha256,
                      const uint32_t* rng, const Measurements& m,
                      const std::vector<ThresholdCheck>& checks) {
    std::string s = "{\"t\":" + json_number(static_cast<double>(t)) + ",\"file\":" +
                    json_string(file) + ",\"sha256\":" + json_string(sha256);
    if (rng) {
        s += ",\"rng\":[";
        for (int q = 0; q < 33; ++q) {
            if (q) s += ",";
            s += std::to_string(rng[q]);
        }
        s += "]";
    }
    s += fields_json(m);
    s += ",\"thresholdChecks\":[";
    for (std::size_t c = 0; c < checks.size(); ++c) {
        if (c) s += ",";
        s += "{\"factor\":" + json_number(checks[c].factor) +
             ",\"threshold\":" + json_number(checks[c].threshold) +
             ",\"mass\":" + json_number(checks[c].mass) +
             ",\"area\":" + json_number(static_cast<double>(checks[c].area)) +
             ",\"density\":" + json_number(checks[c].density) + "}";
    }
    return s + "]}";
}

}  // namespace lgca
