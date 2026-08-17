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

template <class Subdivision>
void feed_subdivision(Subdivision& subdivision, const MulticoloredSimplices& mc_simplices) {
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
    feed_subdivision(subdivision, mc_simplices);
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
    feed_subdivision(subdivision, mc_simplices);
    return mc_simplices;
}

} // namespace delaunay_interfaces::detail
