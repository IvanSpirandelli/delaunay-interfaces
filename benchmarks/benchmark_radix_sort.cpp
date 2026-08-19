// A/B benchmark: the production radix bucket sort (detail::sort_bucket,
// src/filtration_sort.hpp) vs std::sort with the equivalent (value, ids)
// comparator, on real filtration buckets. See benchmarks/RADIX_SORT.md for a
// worked example of the algorithm.
//
// Inputs are the same deterministic point clouds as benchmark_speed.cpp
// (same seeds), plain-Delaunay scenario. The pipeline's finalized buckets are
// shuffled once with a fixed seed to stand in for emission order, and each
// timed repeat sorts a fresh copy of that shuffled vector. The shuffled
// buckets are already deduplicated (the raw emission stream also carries
// exact duplicates, which affect only the sizes, not the comparison), and
// both algorithms are verified byte-identical on every repeat.
//
// Usage: benchmark_radix_sort [--sizes 20000,50000] [--repeats 5]

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <delaunay_interfaces/interface_generation.hpp>
#include "filtration_sort.hpp"
#include "subdivision_driver.hpp"

using namespace delaunay_interfaces;

namespace {

// Seed derived from (n, color mode), identical to benchmark_speed.cpp.
std::pair<Points, ColorLabels> generate_input(size_t n, bool two_colors) {
    std::mt19937 rng(static_cast<unsigned>(1000003 * n + (two_colors ? 1 : 0)));
    std::uniform_real_distribution<double> coord(0.0, 1.0);
    Points points(n);
    ColorLabels colors(n);
    for (size_t i = 0; i < n; ++i) {
        points[i] = Point3D(coord(rng), coord(rng), coord(rng));
        colors[i] = two_colors ? static_cast<int>(i % 2 + 1) : static_cast<int>(i + 1);
    }
    return {std::move(points), std::move(colors)};
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

template <typename Entry, typename FullLess, typename IdLess>
void bench_bucket(const std::string& label, const std::vector<Entry>& sorted,
                  int repeats, FullLess full_less, IdLess id_less) {
    std::vector<Entry> shuffled = sorted;
    std::mt19937 rng(42);
    std::shuffle(shuffled.begin(), shuffled.end(), rng);

    std::vector<double> std_ms, radix_ms;
    for (int rep = 0; rep < repeats; ++rep) {
        std::vector<Entry> a = shuffled;
        auto t0 = std::chrono::steady_clock::now();
        std::sort(a.begin(), a.end(), full_less);
        auto t1 = std::chrono::steady_clock::now();
        std_ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

        std::vector<Entry> b = shuffled;
        t0 = std::chrono::steady_clock::now();
        detail::sort_bucket(b, full_less, id_less);
        t1 = std::chrono::steady_clock::now();
        radix_ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

        if (std::memcmp(a.data(), b.data(), a.size() * sizeof(Entry)) != 0) {
            std::cerr << "MISMATCH: radix output differs from std::sort (" << label << ")\n";
            std::exit(1);
        }
    }

    const double sm = median(std_ms), rm = median(radix_ms);
    std::cout << "  " << std::left << std::setw(22) << label << std::right
              << " | " << std::setw(9) << sorted.size()
              << " | " << std::setw(11) << std::fixed << std::setprecision(2) << sm
              << " | " << std::setw(9) << rm
              << " | " << std::setw(6) << std::setprecision(2) << sm / rm << "x\n";
}

} // namespace

int main(int argc, char** argv) {
    std::vector<size_t> sizes = {20000, 50000};
    int repeats = 5;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--sizes" && i + 1 < argc) {
            sizes.clear();
            std::stringstream ss(argv[++i]);
            std::string tok;
            while (std::getline(ss, tok, ',')) sizes.push_back(std::stoul(tok));
        } else if (arg == "--repeats" && i + 1 < argc) {
            repeats = std::stoi(argv[++i]);
        } else {
            std::cerr << "Usage: benchmark_radix_sort [--sizes n1,n2,...] [--repeats k]\n";
            return 1;
        }
    }

    std::cout << "Radix bucket sort vs std::sort on shuffled pipeline buckets"
              << " (median of " << repeats << " repeats;\nbuckets under "
              << "10000 entries take sort_bucket's std::sort fallback)\n\n"
              << "  scenario / bucket       |   entries | std::sort ms |  radix ms | radix speedup\n"
              << "  ------------------------|-----------|--------------|-----------|--------------\n";

    for (bool two_colors : {true, false}) {
        const std::string mode = two_colors ? "two" : "distinct";
        for (size_t n : sizes) {
            auto [points, colors] = generate_input(n, two_colors);
            auto [verts, flat, mc, vai] = detail::compute_subdivision_flat(
                points, colors, Radii{}, false, false);

            const std::string tag = std::to_string(n) + "/" + mode;
            using D1 = FlatFiltration::Dim1Entry;
            using D2 = FlatFiltration::Dim2Entry;
            using D3 = FlatFiltration::Dim3Entry;

            bench_bucket(tag + " dim1", flat.dim1, repeats,
                [](const D1& a, const D1& b) {
                    if (a.value != b.value) return a.value < b.value;
                    return a.v < b.v;
                },
                [](const D1& a, const D1& b) { return a.v < b.v; });

            bench_bucket(tag + " dim2", flat.dim2, repeats,
                [](const D2& a, const D2& b) {
                    if (a.value != b.value) return a.value < b.value;
                    if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
                    return a.v[1] < b.v[1];
                },
                [](const D2& a, const D2& b) {
                    if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
                    return a.v[1] < b.v[1];
                });

            bench_bucket(tag + " dim3", flat.dim3, repeats,
                [](const D3& a, const D3& b) {
                    if (a.value != b.value) return a.value < b.value;
                    if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
                    if (a.v[1] != b.v[1]) return a.v[1] < b.v[1];
                    return a.v[2] < b.v[2];
                },
                [](const D3& a, const D3& b) {
                    if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
                    if (a.v[1] != b.v[1]) return a.v[1] < b.v[1];
                    return a.v[2] < b.v[2];
                });
        }
    }

    return 0;
}
