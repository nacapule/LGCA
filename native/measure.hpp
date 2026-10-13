// native/measure.hpp — the sweeps' per-sample observables, measured natively.
//
// A port of measurements(occ, W, H, factor), the JavaScript function that the density
// and alpha sweep's runner called at every sample time and whose results fill
// run.json's `series` rows and `snaps`. Given the same lattice it returns the same
// numbers, bit for bit, and writes them as the same JSON text (both are checked
// against the JavaScript function and against the saved runs).
//
// What keeps it bit-identical (docs/FIDELITY.md):
//   - the same arithmetic types as JavaScript: a site's count n is an Int32Array element
//     (it wraps around like a 32-bit integer), every total is a double;
//   - every sum in the same order: sites in index order, channels k = 0..5 within a site,
//     cluster sites in the order the breadth-first search reaches them;
//   - the search visits the 6 neighbours in the same order (NB below) and keeps the
//     first cluster of largest mass (a later cluster replaces it only if heavier);
//   - Math.hypot as V8 computes it (js_hypot below), not std::hypot, which is a
//     different algorithm and can differ in the last bit;
//   - numbers written as JavaScript's JSON.stringify writes them (json_number).
// The histogram and the edge profiles of the JavaScript function are not ported: the
// sweep deletes the histogram before saving, and asks for no profiles.
//
// Lattice layout, as everywhere in native/: occ[(j*W + i)*6 + k] = particles at site
// (i, j) moving in channel k; site index x = j*W + i.

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lgca {

// One sample's observables, the fields of a run.json `series` row (same names, same
// order). With n = particles at a site, p = its momentum sum_k n_k c_k (c_k the unit
// velocity of channel k), A = W*H sites:
struct Measurements {
    double N = 0;                   // particles, sum of n
    double mean = 0;                // N / A
    int32_t nmax = 0;               // largest n
    int32_t kmax = 0;               // largest single-channel count
    double siteFraction = 0;        // nmax / N             (0 when N = 0)
    double channelFraction = 0;     // kmax / N             (0 when N = 0)
    double effectiveArea = 0;       // N^2 / sum n^2        (0 when every n = 0)
    int64_t occupied = 0;           // sites with n != 0
    double polar = 0;               // |sum of p| / N       (0 when N = 0)
    double localPolar = 0;          // mean of |p| / n over occupied sites
    double weightedLocalPolar = 0;  // sum of |p| / N       (0 when N = 0)
    double threshold = 0;           // cluster threshold max(1, ceil(factor * N / A))
    int64_t clusters = 0;           // clusters: connected sets of sites with n >= threshold
    double clusterMass = 0;         // particles in the heaviest cluster
    int64_t clusterArea = 0;        // its sites
    double clusterDensity = 0;      // clusterMass / clusterArea   (0 without clusters)
    double clusterFraction = 0;     // clusterMass / N             (0 when N = 0)
    double clusterPolar = 0;        // |sum of its p| / clusterMass (0 without clusters)
};

// One entry of a snap's `thresholdChecks`: the heaviest cluster when the threshold uses
// another factor (the sweep records factors 1.5 and 4).
struct ThresholdCheck {
    double factor;
    double threshold;
    double mass;        // its clusterMass
    int64_t area;       // its clusterArea
    double density;     // its clusterDensity
};

// The factors of a snap's thresholdChecks, in the sweep's order.
constexpr double THRESHOLD_CHECK_FACTORS[2] = {1.5, 4};

// The island (defined for the island study on 2026-09-30): the site holding the most
// particles, together with every occupied site (n >= 1) connected to it through occupied
// neighbours, using the same 6 periodic neighbours as the clusters. No density threshold:
// any particle on a site links that site. The busiest site is the first one in index
// order when several share the largest n. An empty lattice has no island (seed -1,
// everything 0 apart from densest2Mass, NaN).
//
// The islands of a lattice are all such connected sets of occupied sites; `components`
// counts them. Three of them are recorded besides the island, each by its particles
// (mass) and sites (area):
//   - the heaviest: the most particles (largestMass, largestArea);
//   - the densest: the most particles per site, mass / area (densestMass, densestArea);
//     often a single site, a pile with no occupied neighbour;
//   - the densest of at least two sites (densest2Mass, densest2Area); NaN and 0 when
//     every island is a single site.
// When islands tie, the one whose lowest site index is smallest is chosen: the order
// and rule of results/islands-2026-10-06/build/heaviest.py, which lists the islands by
// their lowest site index and keeps the first of equal values (numpy's argmax). Densities
// are compared as the doubles mass / area, as there, so two islands whose quotients
// round to the same double tie. Masses are double sums, exact while an island holds
// fewer than 2^53 particles; above that, this search and heaviest.py round differently.
// The island is usually the heaviest, but need not be: a tall pile with little around it
// can hold fewer particles than a wide, lower island.
//
// Positions use unwrapped coordinates: the search follows the island across the
// periodic edges and gives each site the coordinates of the path that reached it. If an
// island reaches one of its own sites along two paths whose coordinates differ by a
// whole box, it closes around the box (`wraps`): it spans the lattice, and its centre
// and radius of gyration are undefined (NaN).
struct Island {
    int32_t seed = -1;          // busiest site's index j*W + i; -1 when N = 0
    double mass = 0;            // particles in the island
    int64_t area = 0;           // its sites
    double polar = 0;           // |sum of its p| / mass   (0 when mass = 0)
    bool wraps = false;         // it closes around the periodic box
    double ci = 0, cj = 0;      // mass-weighted centre in lattice coordinates (i, j), in [0, W) x [0, H)
    double rg = 0;              // mass-weighted radius of gyration, in lattice spacings
    int64_t components = 0;     // number of islands in the lattice
    double largestMass = 0;     // the heaviest island's mass
    int64_t largestArea = 0;    //   and sites
    double densestMass = 0;     // the densest island's mass
    int64_t densestArea = 0;    //   and sites
    double densest2Mass = std::numeric_limits<double>::quiet_NaN();   // the densest island
    int64_t densest2Area = 0;   //   of at least two sites: mass and sites (NaN, 0: none)
};

// densest2Area as a number to write: NaN (JSON null, like densest2Mass) when there is no
// island of at least two sites.
inline double island_densest2_area(const Island& s) {
    return s.densest2Area > 0 ? static_cast<double>(s.densest2Area)
                              : std::numeric_limits<double>::quiet_NaN();
}

// Measures lattices of one size. Its work arrays are made once, at their full size, so
// measure() allocates no memory.
class Measurer {
public:
    Measurer(int W, int H);

    int width() const { return W_; }
    int height() const { return H_; }

    // The JavaScript measurements(occ, W, H, factor), without histogram and profiles.
    // occ holds W*H*6 values. The sweeps' series use factor 2.
    Measurements measure(const int32_t* occ, double factor = 2);

    // The sweep's thresholdChecks: measure(occ, f) for f = 1.5 and 4, reduced to
    // {factor, threshold, mass, area, density}.
    std::vector<ThresholdCheck> threshold_checks(const int32_t* occ);

    // The island of the lattice and the count of all islands (see Island). Independent
    // of measure(): it reads occ itself and changes nothing measure() returns. Throws
    // std::invalid_argument if a channel count is negative, std::overflow_error if a
    // site's total does not fit in 32 bits.
    Island island(const int32_t* occ);

private:
    int W_, H_;
    std::vector<int32_t> n_;       // particles per site (an Int32Array in JavaScript)
    std::vector<double> px_, py_;  // momentum per site
    std::vector<uint8_t> seen_;    // the cluster search's visited marks
    std::vector<int32_t> queue_;   // the cluster search's queue
    std::vector<int32_t> best_;    // sites of the heaviest cluster, in queue order
    std::vector<int32_t> ui_, uj_; // the island search's unwrapped coordinates
};

// An Island as JSON fields: ,"seed":...,"islandMass":...,... (each preceded by a comma).
std::string island_fields_json(const Island& s);

// Math.hypot(a, b) as V8 (Node, Chrome) computes it: both values divided by the larger
// magnitude, the squares added with Kahan's compensated summation, then the square root
// times that magnitude. Basic IEEE operations only, so the result is the same double on
// any machine (compiled with -ffp-contract=off).
double js_hypot(double a, double b);

// A number as JavaScript's JSON.stringify writes it: the shortest decimal that reads back
// as the same double, in JavaScript's layout (1228, 0.15160493827160493, 1e-7, 1e+21);
// -0 as 0; NaN and infinities as null.
std::string json_number(double v);

// A string as a JSON string literal, quotes included, escaped as JSON.stringify escapes
// it ("state-t100.bin.gz" stays as it is; a quote, a backslash or a control character
// gets its escape).
std::string json_string(const std::string& s);

// A run.json `series` row: {"t":t,"N":...,...,"clusterPolar":...} on one line, the text
// JSON.stringify({t, ...measurements}) gives.
std::string series_json(long t, const Measurements& m);

// A run.json `snaps` entry: {"t","file","sha256","rng",<the series fields>,
// "thresholdChecks":[...]}, the text JSON.stringify gives for the sweep's snap object.
// rng points to the 33 generator words (32 state words, then the index), or is null
// when the lattice came without them: the "rng" field is then left out.
std::string snap_json(long t, const std::string& file, const std::string& sha256,
                      const uint32_t* rng, const Measurements& m,
                      const std::vector<ThresholdCheck>& checks);

}  // namespace lgca
