#pragma once

#include "delaunay_interfaces/types.hpp"
#include <cmath>

namespace test_helpers {

// Single tolerance rule for all floating-point comparisons: relative for
// large magnitudes with an absolute floor near zero.
inline bool approx_equal(double a, double b, double rel = 1e-9, double abs_tol = 1e-12) {
    return std::abs(a - b) <= std::max(abs_tol, rel * std::max(std::abs(a), std::abs(b)));
}

struct SimplexCounts {
    size_t vertices = 0;
    size_t edges = 0;
    size_t triangles = 0;

    size_t total() const { return vertices + edges + triangles; }
    long euler_characteristic() const {
        return static_cast<long>(vertices) - static_cast<long>(edges) + static_cast<long>(triangles);
    }
};

inline SimplexCounts count_simplices(const delaunay_interfaces::Filtration& filtration) {
    SimplexCounts counts;
    for (const auto& [simplex, value] : filtration) {
        switch (simplex.size()) {
        case 1: ++counts.vertices; break;
        case 2: ++counts.edges; break;
        case 3: ++counts.triangles; break;
        default: break;
        }
    }
    return counts;
}

} // namespace test_helpers
