// Benchmark suite: random point clouds in the unit cube at several densities,
// two color modes (all-distinct: subdivision-heavy; 2 colors: construction-
// heavy). Per point cloud the full unweighted Delaunay interface is computed
// once and its median filtration value L taken as the data-derived length
// scale (the maximum is dominated by boundary slivers and yields near-full
// alpha complexes); alpha complexes are then benchmarked at radius L/2 and
// L/4 — uniform (scalar radius) and weighted (heterogeneous radii drawn
// around the scale, so hidden points and weighted predicates are exercised).
//
// Per scenario the collect phase (triangulation + multicolored simplex
// collection) and the full pipeline are timed separately (median of repeats),
// and output invariants (counts + filtration-value sum) are recorded so
// optimization commits can prove output-identical results.
//
// Usage: benchmark_speed [--sizes 1000,5000] [--repeats 3] [--out results.json]
// Human-readable table goes to stdout; JSON records to --out when given.

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <delaunay_interfaces/interface_generation.hpp>

using namespace delaunay_interfaces;

namespace {

struct InputData {
    Points points;
    ColorLabels colors;
};

// Seed derived from (n, color mode) so every branch/commit benchmarks
// byte-identical inputs.
InputData generate_input(size_t n, bool two_colors) {
    std::mt19937 rng(static_cast<unsigned>(1000003 * n + (two_colors ? 1 : 0)));
    std::uniform_real_distribution<double> coord(0.0, 1.0);

    InputData data;
    data.points.resize(n);
    data.colors.resize(n);
    for (size_t i = 0; i < n; ++i) {
        data.points[i] = Point3D(coord(rng), coord(rng), coord(rng));
        data.colors[i] = two_colors ? static_cast<int>(i % 2 + 1) : static_cast<int>(i + 1);
    }
    return data;
}

// Heterogeneous radii around scale s: U(0.8 s, 1.2 s), fixed seed.
Radii generate_radii(size_t n, double scale, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.8 * scale, 1.2 * scale);
    Radii radii(n);
    for (auto& r : radii) r = dist(rng);
    return radii;
}

struct Invariants {
    size_t n_vertices = 0;
    size_t n_simplices = 0;
    size_t n_tets = 0;
    size_t n_free_triangles = 0;
    size_t n_free_edges = 0;
    double filtration_sum = 0.0;
    double median_filtration = 0.0;
};

struct Result {
    double collect_ms = 0.0;
    double total_ms = 0.0;
    Invariants inv;
};

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// One benchmarked scenario: `collect` and `total` run the two measured
// phases on identical inputs; the pipeline result feeds the invariants.
template <class CollectFn, class TotalFn>
Result run_scenario(int repeats, CollectFn collect, TotalFn total) {
    Result res;
    std::vector<double> collect_times, total_times;

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        auto simplices = collect();
        auto t1 = std::chrono::steady_clock::now();
        collect_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

        t0 = std::chrono::steady_clock::now();
        auto [verts, filtration, mc, vai] = total();
        t1 = std::chrono::steady_clock::now();
        total_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

        if (rep == 0) {
            res.inv.n_vertices = verts.size();
            res.inv.n_simplices = filtration.size();
            res.inv.n_tets = mc.tetrahedra.size();
            res.inv.n_free_triangles = mc.free_triangles.size();
            res.inv.n_free_edges = mc.free_edges.size();
            std::vector<double> values;
            values.reserve(filtration.size());
            for (const auto& [s, val] : filtration) {
                res.inv.filtration_sum += val;
                values.push_back(val);
            }
            if (!values.empty()) res.inv.median_filtration = median(std::move(values));
        }
    }

    res.collect_ms = median(collect_times);
    res.total_ms = median(total_times);
    return res;
}

void print_row(size_t n, const std::string& mode, const std::string& scenario, const Result& r) {
    std::cout << "  " << std::setw(7) << n
              << " | " << std::left << std::setw(8) << mode
              << " | " << std::setw(25) << scenario << std::right
              << " | " << std::setw(9) << std::fixed << std::setprecision(1) << r.collect_ms
              << " | " << std::setw(9) << r.total_ms
              << " | " << std::setw(8) << r.inv.n_vertices
              << " | " << std::setw(9) << r.inv.n_simplices << "\n";
}

void append_json(std::ostringstream& json, bool& first, size_t n, const std::string& mode,
                 const std::string& scenario, double L, const Result& r) {
    if (!first) json << ",\n";
    first = false;
    json << "    {\"n\": " << n
         << ", \"colors\": \"" << mode << "\""
         << ", \"scenario\": \"" << scenario << "\""
         << ", \"L\": " << std::setprecision(17) << L
         << ", \"collect_ms\": " << std::setprecision(6) << r.collect_ms
         << ", \"total_ms\": " << r.total_ms
         << ", \"n_vertices\": " << r.inv.n_vertices
         << ", \"n_simplices\": " << r.inv.n_simplices
         << ", \"n_tets\": " << r.inv.n_tets
         << ", \"n_free_triangles\": " << r.inv.n_free_triangles
         << ", \"n_free_edges\": " << r.inv.n_free_edges
         << ", \"filtration_sum\": " << std::setprecision(17) << r.inv.filtration_sum
         << "}";
}

} // namespace

int main(int argc, char** argv) {
    std::vector<size_t> sizes = {1000, 5000, 20000, 50000};
    int repeats = 3;
    std::string out_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--sizes" && i + 1 < argc) {
            sizes.clear();
            std::stringstream ss(argv[++i]);
            std::string tok;
            while (std::getline(ss, tok, ',')) sizes.push_back(std::stoul(tok));
        } else if (arg == "--repeats" && i + 1 < argc) {
            repeats = std::stoi(argv[++i]);
        } else if (arg == "--out" && i + 1 < argc) {
            out_path = argv[++i];
        } else {
            std::cerr << "Usage: benchmark_speed [--sizes n1,n2,...] [--repeats k] [--out file.json]\n";
            return 1;
        }
    }

    InterfaceGenerator gen;
    std::ostringstream json;
    bool first = true;
    json << "[\n";

    std::cout << "Benchmark: unit-cube point clouds (median of " << repeats << " repeats)\n\n"
              << "  N pts   | Colors   | Scenario                  | col'ct ms |  total ms | Vertices | Simplices\n"
              << "  --------|----------|---------------------------|-----------|-----------|----------|----------\n";

    for (bool two_colors : {false, true}) {
        const std::string mode = two_colors ? "two" : "distinct";
        for (size_t n : sizes) {
            auto data = generate_input(n, two_colors);
            const auto& pts = data.points;
            const auto& cols = data.colors;

            // Unweighted Delaunay; its median filtration value is the length scale L.
            auto res = run_scenario(repeats,
                [&] { return gen.collect_multicolored_simplices(pts, cols); },
                [&] { return compute_barycentric_subdivision_and_filtration(pts, cols, Radii{}, false); });
            const double L = res.inv.median_filtration;
            print_row(n, mode, "delaunay", res);
            append_json(json, first, n, mode, "delaunay", L, res);

            const unsigned radii_seed = static_cast<unsigned>(2000003 * n + (two_colors ? 1 : 0));
            for (double frac : {0.5, 0.25}) {
                const std::string suffix = frac == 0.5 ? "half" : "quarter";
                const double scale = L * frac;
                auto radii = generate_radii(n, scale, radii_seed + (frac == 0.5 ? 0 : 1));

                res = run_scenario(repeats,
                    [&] { return gen.collect_multicolored_simplices(pts, cols, radii, false); },
                    [&] { return compute_barycentric_subdivision_and_filtration(pts, cols, radii, false); });
                print_row(n, mode, "delaunay_weighted_" + suffix, res);
                append_json(json, first, n, mode, "delaunay_weighted_" + suffix, L, res);

                res = run_scenario(repeats,
                    [&] { return gen.collect_multicolored_simplices(pts, cols, radii, true); },
                    [&] { return compute_barycentric_subdivision_and_filtration(pts, cols, radii, true); });
                print_row(n, mode, "alpha_weighted_" + suffix, res);
                append_json(json, first, n, mode, "alpha_weighted_" + suffix, L, res);

                res = run_scenario(repeats,
                    [&] { return gen.collect_multicolored_simplices(pts, cols, scale); },
                    [&] { return compute_barycentric_subdivision_and_filtration(pts, cols, scale); });
                print_row(n, mode, "alpha_uniform_" + suffix, res);
                append_json(json, first, n, mode, "alpha_uniform_" + suffix, L, res);
            }
            std::cout << "  --------|----------|---------------------------|-----------|-----------|----------|----------\n";
        }
    }

    json << "\n]\n";

    if (!out_path.empty()) {
        std::ofstream out(out_path);
        if (!out) {
            std::cerr << "Cannot write " << out_path << "\n";
            return 1;
        }
        out << json.str();
        std::cout << "\nJSON written to " << out_path << "\n";
    }

    return 0;
}
