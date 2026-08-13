// Benchmark: random point clouds with all distinct colors.
// Measures wall-clock time for compute_barycentric_subdivision_and_filtration.

#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>
#include <delaunay_interfaces/interface_generation.hpp>
#include <delaunay_interfaces/barycentric_subdivision.hpp>

using namespace delaunay_interfaces;

struct TimingResult {
    size_t n_vertices;
    size_t n_simplices;
    double time_ms;
};

struct InputData {
    Points points;
    ColorLabels colors;
    Radii radii;
};

InputData generate_input(size_t n, unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-10.0, 10.0);
    std::uniform_real_distribution<double> radius(0.5, 2.0);

    InputData data;
    data.points.resize(n);
    data.colors.resize(n);
    data.radii.resize(n);

    for (size_t i = 0; i < n; ++i) {
        data.points[i] = Point3D(coord(rng), coord(rng), coord(rng));
        data.colors[i] = static_cast<int>(i + 1);
        data.radii[i] = radius(rng);
    }
    return data;
}

TimingResult run_test(const InputData& data, bool weighted, bool alpha) {
    auto start = std::chrono::high_resolution_clock::now();

    auto [verts, filtration, mc_simplices, vai] = compute_barycentric_subdivision_and_filtration(
        data.points, data.colors, data.radii, weighted, alpha);

    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();

    return {verts.size(), filtration.size(), ms};
}

void print_row(const std::string& n, const std::string& mode, const TimingResult& r) {
    std::cout << "  " << std::setw(8) << n
              << " | " << std::left << std::setw(12) << mode << std::right
              << " | " << std::setw(8) << r.n_vertices
              << " | " << std::setw(10) << r.n_simplices
              << " | " << std::setw(8) << std::fixed << std::setprecision(1) << r.time_ms
              << " ms\n";
}

int main() {
    std::cout << "Benchmark: Random Point Clouds (all distinct colors)\n\n";

    std::cout << "  N points |     Mode     | Vertices |  Simplices |      Time\n";
    std::cout << "  ---------|--------------|----------|------------|----------\n";

    for (size_t n : {10, 100, 1000, 10000, 25000, 50000}) {
        auto data = generate_input(n);

        print_row(std::to_string(n), "Delaunay", run_test(data, true, false));
        print_row("", "Alpha", run_test(data, true, true));

        std::cout << "  ---------|--------------|----------|------------|----------\n";
    }

    return 0;
}
