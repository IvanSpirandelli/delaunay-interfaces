// Measurement driver for the README performance figure: four-color random
// clouds, four pipeline flavors, N repetitions each, all raw times emitted
// as JSON so the plot can show mean + run-to-run spread.
//
// Build (from repo root):
//   c++ -std=c++17 -O3 -DNDEBUG -Iinclude -I/opt/homebrew/include \
//       -I/opt/homebrew/include/eigen3 benchmarks/visualization/readme_bench.cpp \
//       build/libdelaunay_interfaces.a -L/opt/homebrew/lib -lgmp -lmpfr \
//       -o benchmarks/visualization/readme_bench
// Run:  readme_bench <reps> <out.json>
#include <delaunay_interfaces/interface_generation.hpp>

#include <chrono>
#include <cstdio>
#include <functional>
#include <random>
#include <vector>

using namespace delaunay_interfaces;

namespace {

struct SizeSpec {
    size_t n;
    double L; // length scale: median filtration of the plain Delaunay run
};
// L values from benchmarks/results/restructure-examples_c71964c.json.
const std::vector<SizeSpec> kSizes = {
    {1000, 0.09431918495730701},
    {5000, 0.0548053427061956},
    {20000, 0.03444298001718487},
    {50000, 0.025361934564134252},
};
constexpr int kColors = 4;

double run_ms(const std::function<InterfaceSurface()>& f, size_t& sink) {
    auto t0 = std::chrono::steady_clock::now();
    InterfaceSurface s = f();
    auto t1 = std::chrono::steady_clock::now();
    sink += s.filtration.size() + s.vertices.size();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: readme_bench <reps> <out.json>\n");
        return 1;
    }
    const int reps = std::atoi(argv[1]);
    std::FILE* out = std::fopen(argv[2], "w");
    if (!out) return 1;

    InterfaceGenerator gen;
    size_t sink = 0;
    std::fprintf(out, "{\n  \"colors\": %d,\n  \"reps\": %d,\n  \"groups\": {\n", kColors, reps);

    const char* group_names[4] = {
        "full Delaunay, unweighted", "full Delaunay, weighted",
        "alpha complex, weighted", "alpha complex, uniform radius"};
    bool first_group = true;
    for (int g = 0; g < 4; ++g) {
        if (!first_group) std::fprintf(out, ",\n");
        first_group = false;
        std::fprintf(out, "    \"%s\": {", group_names[g]);
        bool first_size = true;
        for (const auto& spec : kSizes) {
            const size_t n = spec.n;
            std::mt19937 rng(static_cast<unsigned>(1000003 * n + 4));
            std::uniform_real_distribution<double> coord(0.0, 1.0);
            Points pts(n);
            ColorLabels cols(n);
            for (size_t i = 0; i < n; ++i) {
                pts[i] = Point3D(coord(rng), coord(rng), coord(rng));
                cols[i] = static_cast<int>(i % kColors + 1);
            }
            const double scale = spec.L / 2;
            Radii radii(n);
            std::mt19937 rrng(static_cast<unsigned>(2000003 * n + 4));
            std::uniform_real_distribution<double> rdist(0.8 * scale, 1.2 * scale);
            for (auto& r : radii) r = rdist(rrng);

            if (!first_size) std::fprintf(out, ",");
            first_size = false;
            std::fprintf(out, "\n      \"%zu\": [", n);
            for (int rep = 0; rep < reps; ++rep) {
                double ms = 0.0;
                switch (g) {
                    case 0: ms = run_ms([&] { return gen.compute_interface_surface(pts, cols); }, sink); break;
                    case 1: ms = run_ms([&] { return gen.compute_interface_surface(pts, cols, radii, false); }, sink); break;
                    case 2: ms = run_ms([&] { return gen.compute_interface_surface(pts, cols, radii, true); }, sink); break;
                    case 3: ms = run_ms([&] { return gen.compute_interface_surface(pts, cols, scale); }, sink); break;
                }
                std::fprintf(out, "%s%.3f", rep ? ", " : "", ms);
            }
            std::fprintf(out, "]");
            std::fprintf(stderr, "%s n=%zu done\n", group_names[g], n);
        }
        std::fprintf(out, "\n    }");
    }
    std::fprintf(out, "\n  }\n}\n");
    std::fclose(out);
    std::fprintf(stderr, "sink %zu\n", sink);
    return 0;
}
