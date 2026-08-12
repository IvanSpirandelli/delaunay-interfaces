#pragma once

#include "types.hpp"

namespace delaunay_interfaces {

class InterfaceGenerator {
public:
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        bool weighted = true,
        bool alpha = true,
        bool lower_star = false
    ) const;

    // Radii-free overload: plain (unweighted) Delaunay complex.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    [[nodiscard]] Tetrahedra get_multicolored_tetrahedra(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        bool weighted = true,
        bool alpha = true
    ) const;

    [[nodiscard]] Tetrahedra get_multicolored_tetrahedra(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    // All multicolored simplices of the weighted alpha complex:
    // tetrahedra + free triangles + free edges not covered by tetrahedra.
    [[nodiscard]] MulticoloredSimplices get_multicolored_simplices_weighted_alpha(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    ) const;

private:
    Tetrahedra get_multicolored_tetrahedra_delaunay(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    Tetrahedra get_multicolored_tetrahedra_weighted_delaunay(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    ) const;

    Tetrahedra get_multicolored_tetrahedra_weighted_alpha(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    ) const;

    bool is_multicolored(const Tetrahedron& tet, const ColorLabels& color_labels) const;
};

[[nodiscard]] std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> get_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted = true,
    bool alpha = true,
    bool lower_star = false
);

// Simplified 2-color surface: each vertex = midpoint of one cross-color atom pair.
// Requires exactly 2 distinct colors in the input.
[[nodiscard]] SimplifiedSurface get_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted = true,
    bool alpha = true
);

} // namespace delaunay_interfaces
