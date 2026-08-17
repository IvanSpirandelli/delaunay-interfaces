// Internal header: shared driver for feeding multicolored simplices of the
// chosen complex into a subdivision object (BarycentricSubdivision or
// SimplifiedSubdivision, which share the process_* interface).
#pragma once

#include "delaunay_interfaces/interface_generation.hpp"
#include <stdexcept>

namespace delaunay_interfaces::detail {

inline void validate_inputs(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) {
    if (points.size() != color_labels.size()) {
        throw std::invalid_argument("Each point must have a corresponding color_label");
    }
    // Non-empty radii select the weighted complex and must cover every point.
    if (!radii.empty() && radii.size() != points.size()) {
        throw std::invalid_argument("Each point must have an assigned radius for weighted complexes");
    }
}

// Exact per-simplex vertex (= multicolored face; an upper bound on new
// vertices, since adjacent simplices share them) and non-singleton
// filtration-entry counts of the barycentric subdivision, per chromatic
// partition shape. Derived by enumerating the face poset of each shape and
// verified against process_simplex output:
//   tet 2-2: 9/24, 3-1: 7/18, 2-1-1: 10/29, 1-1-1-1: 11/34,
//   triangle 2-1: 3/2, 1-1-1: 4/3, edge 1-1: 1/0.
struct SubdivisionCapacity {
    size_t vertices = 0;
    size_t entries = 0;
};

inline SubdivisionCapacity subdivision_capacity(
    const ColorLabels& color_labels,
    const MulticoloredSimplices& mc_simplices
) {
    SubdivisionCapacity cap;
    for (const auto& tet : mc_simplices.tetrahedra) {
        // Distinct colors and largest same-color group among the 4 labels,
        // via fixed arrays (no allocation).
        int colors[4];
        int counts[4];
        int k = 0;
        for (int v : tet) {
            int c = color_labels[v];
            int j = 0;
            while (j < k && colors[j] != c) ++j;
            if (j == k) {
                colors[k] = c;
                counts[k] = 0;
                ++k;
            }
            ++counts[j];
        }
        if (k == 2) {
            const bool two_two = counts[0] == 2;
            cap.vertices += two_two ? 9 : 7;
            cap.entries += two_two ? 24 : 18;
        } else if (k == 3) {
            cap.vertices += 10;
            cap.entries += 29;
        } else if (k == 4) {
            cap.vertices += 11;
            cap.entries += 34;
        }
    }
    for (const auto& tri : mc_simplices.free_triangles) {
        const bool two_colors = color_labels[tri[0]] == color_labels[tri[1]]
            || color_labels[tri[0]] == color_labels[tri[2]]
            || color_labels[tri[1]] == color_labels[tri[2]];
        cap.vertices += two_colors ? 3 : 4;
        cap.entries += two_colors ? 2 : 3;
    }
    cap.vertices += mc_simplices.free_edges.size();
    return cap;
}

template <class Subdivision>
void feed_subdivision(
    Subdivision& subdivision,
    const ColorLabels& color_labels,
    const MulticoloredSimplices& mc_simplices
) {
    auto cap = subdivision_capacity(color_labels, mc_simplices);
    subdivision.reserve(cap.vertices, cap.entries);
    for (const auto& tet : mc_simplices.tetrahedra) {
        subdivision.process_tetrahedron(tet);
    }
    for (const auto& tri : mc_simplices.free_triangles) {
        subdivision.process_simplex({tri.begin(), tri.end()});
    }
    for (const auto& edge : mc_simplices.free_edges) {
        subdivision.process_simplex({edge.begin(), edge.end()});
    }
}

template <class Subdivision>
MulticoloredSimplices run_subdivision(
    Subdivision& subdivision,
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool alpha
) {
    InterfaceGenerator generator;
    MulticoloredSimplices mc_simplices = generator.collect_multicolored_simplices(
        points, color_labels, radii, alpha);
    feed_subdivision(subdivision, color_labels, mc_simplices);
    return mc_simplices;
}

template <class Subdivision>
MulticoloredSimplices run_subdivision(
    Subdivision& subdivision,
    const Points& points,
    const ColorLabels& color_labels,
    double radius
) {
    InterfaceGenerator generator;
    MulticoloredSimplices mc_simplices = generator.collect_multicolored_simplices(
        points, color_labels, radius);
    feed_subdivision(subdivision, color_labels, mc_simplices);
    return mc_simplices;
}

} // namespace delaunay_interfaces::detail
