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
    const Radii& radii,
    bool weighted
) {
    if (points.size() != color_labels.size()) {
        throw std::invalid_argument("Each point must have a corresponding color_label");
    }
    if (weighted && radii.size() != points.size()) {
        throw std::invalid_argument("Each point must have an assigned radius for weighted complexes");
    }
}

template <class Subdivision>
MulticoloredSimplices run_subdivision(
    Subdivision& subdivision,
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) {
    InterfaceGenerator generator;
    MulticoloredSimplices mc_simplices;

    if (weighted && alpha) {
        mc_simplices = generator.get_multicolored_simplices_weighted_alpha(
            points, color_labels, radii);

        for (const auto& tet : mc_simplices.generating_tetrahedra) {
            subdivision.process_tetrahedron(tet);
        }
        for (const auto& tri : mc_simplices.generating_free_triangles) {
            subdivision.process_simplex(tri);
        }
        for (const auto& edge : mc_simplices.generating_free_edges) {
            subdivision.process_simplex(edge);
        }
    } else {
        auto tetrahedra = generator.get_multicolored_tetrahedra(
            points, color_labels, radii, weighted, alpha);

        for (const auto& tet : tetrahedra) {
            subdivision.process_tetrahedron(tet);
        }
        mc_simplices.generating_tetrahedra = std::move(tetrahedra);
    }

    return mc_simplices;
}

} // namespace delaunay_interfaces::detail
