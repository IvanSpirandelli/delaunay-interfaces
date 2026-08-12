// Shared JSON ground-truth harness. Compares a computed InterfaceSurface
// against reference output produced by the original Julia implementation
// (DelaunayInterfaces.jl/generate_ground_truth.jl), including the full
// simplex structure, not just counts and value multisets.
#pragma once

#include "delaunay_interfaces/interface_generation.hpp"
#include "test_helpers.hpp"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace ground_truth {

using json = nlohmann::json;
using namespace delaunay_interfaces;

struct CaseResult {
    bool passed = false;
    std::string message;
    size_t computed_vertices = 0;
    size_t expected_vertices = 0;
    size_t computed_simplices = 0;
    size_t expected_simplices = 0;
};

inline CaseResult run_case(const json& test_data) {
    CaseResult result;

    const auto& input = test_data["input"];
    const auto& output = test_data["output"];

    Points points;
    for (const auto& p : input["points"]) {
        points.push_back({p[0].get<double>(), p[1].get<double>(), p[2].get<double>()});
    }
    ColorLabels colors;
    for (const auto& c : input["color_labels"]) {
        colors.push_back(c.get<int>());
    }
    Radii radii;
    for (const auto& r : input["radii"]) {
        radii.push_back(r.get<double>());
    }
    bool weighted = input["weighted"].get<bool>();
    bool alpha = input["alpha"].get<bool>();

    result.expected_vertices = output["num_vertices"].get<size_t>();
    result.expected_simplices = output["num_simplices"].get<size_t>();

    InterfaceGenerator generator;
    InterfaceSurface surface;
    try {
        surface = generator.compute_interface_surface(points, colors, radii, weighted, alpha);
    } catch (const std::exception& e) {
        result.message = std::string("Exception: ") + e.what();
        return result;
    }

    result.computed_vertices = surface.vertices.size();
    result.computed_simplices = surface.filtration.size();

    if (result.computed_vertices != result.expected_vertices) {
        result.message = "Vertex count mismatch: " + std::to_string(result.computed_vertices) +
                         " vs expected " + std::to_string(result.expected_vertices);
        return result;
    }
    if (result.computed_simplices != result.expected_simplices) {
        result.message = "Simplex count mismatch: " + std::to_string(result.computed_simplices) +
                         " vs expected " + std::to_string(result.expected_simplices);
        return result;
    }
    if (result.expected_vertices == 0) {
        result.passed = true;
        result.message = "OK";
        return result;
    }

    // Vertex IDs are assigned in creation order, which the C++ port need not
    // share with the Julia reference. Pair vertices by sorted coordinates and
    // derive the index mapping from that pairing.
    using IndexedPoint = std::pair<std::array<double, 3>, size_t>;
    std::vector<IndexedPoint> computed_pts, expected_pts;
    for (size_t i = 0; i < surface.vertices.size(); ++i) {
        const auto& v = surface.vertices[i];
        computed_pts.push_back({{v.x(), v.y(), v.z()}, i});
    }
    {
        size_t i = 0;
        for (const auto& v : output["vertices"]) {
            expected_pts.push_back({{v[0].get<double>(), v[1].get<double>(), v[2].get<double>()}, i++});
        }
    }
    std::sort(computed_pts.begin(), computed_pts.end());
    std::sort(expected_pts.begin(), expected_pts.end());

    std::vector<size_t> computed_to_expected(computed_pts.size());
    for (size_t i = 0; i < computed_pts.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (!test_helpers::approx_equal(computed_pts[i].first[j], expected_pts[i].first[j])) {
                result.message = "Vertex " + std::to_string(i) + " coord " + std::to_string(j) +
                                 " mismatch: " + std::to_string(computed_pts[i].first[j]) +
                                 " vs " + std::to_string(expected_pts[i].first[j]);
                return result;
            }
        }
        computed_to_expected[computed_pts[i].second] = expected_pts[i].second;
    }

    // Compare the filtration as sorted (simplex, value) pairs in the expected
    // vertex numbering. JSON simplices are 1-based (Julia); convert to 0-based.
    using Entry = std::pair<std::vector<int>, double>;
    std::vector<Entry> computed_entries, expected_entries;

    for (const auto& [simplex, val] : surface.filtration) {
        std::vector<int> mapped;
        for (auto id : simplex) {
            mapped.push_back(static_cast<int>(computed_to_expected[static_cast<size_t>(id)]));
        }
        std::sort(mapped.begin(), mapped.end());
        computed_entries.push_back({std::move(mapped), val});
    }
    for (const auto& f : output["filtration"]) {
        std::vector<int> simplex;
        for (const auto& id : f["simplex"]) {
            simplex.push_back(id.get<int>() - 1);
        }
        std::sort(simplex.begin(), simplex.end());
        expected_entries.push_back({std::move(simplex), f["value"].get<double>()});
    }

    auto by_simplex = [](const Entry& a, const Entry& b) { return a.first < b.first; };
    std::sort(computed_entries.begin(), computed_entries.end(), by_simplex);
    std::sort(expected_entries.begin(), expected_entries.end(), by_simplex);

    for (size_t i = 0; i < computed_entries.size(); ++i) {
        if (computed_entries[i].first != expected_entries[i].first) {
            result.message = "Simplex structure mismatch at sorted position " + std::to_string(i);
            return result;
        }
        if (!test_helpers::approx_equal(computed_entries[i].second, expected_entries[i].second)) {
            result.message = "Filtration value mismatch for simplex at sorted position " +
                             std::to_string(i) + ": " + std::to_string(computed_entries[i].second) +
                             " vs " + std::to_string(expected_entries[i].second);
            return result;
        }
    }

    result.passed = true;
    result.message = "OK";
    return result;
}

inline bool report(const std::string& name, const CaseResult& result) {
    std::cout << "Test: " << name << "\n"
              << "  Computed: " << result.computed_vertices << " vertices, "
              << result.computed_simplices << " simplices\n"
              << "  Expected: " << result.expected_vertices << " vertices, "
              << result.expected_simplices << " simplices\n"
              << (result.passed ? "  PASS" : "  FAIL: " + result.message) << "\n";
    return result.passed;
}

// Runs every case in a {name: {input, output}} JSON file. A missing or
// unreadable file is a failure, not a skip.
inline bool run_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "FAIL: Could not open " << filepath << "\n";
        return false;
    }

    json data = json::parse(file);
    bool all_passed = true;
    for (auto& [name, test_data] : data.items()) {
        if (!report(name, run_case(test_data))) {
            all_passed = false;
        }
    }
    return all_passed;
}

// Runs a file whose root is a single {input, output} case.
inline bool run_single_case_file(const std::string& name, const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "FAIL: Could not open " << filepath << "\n";
        return false;
    }
    return report(name, run_case(json::parse(file)));
}

} // namespace ground_truth
