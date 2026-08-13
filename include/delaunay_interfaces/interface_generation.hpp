#pragma once

#include "types.hpp"

namespace delaunay_interfaces {

class InterfaceGenerator {
public:
    // Weightedness is determined by the radii: empty radii build the plain
    // Delaunay complex, non-empty radii the weighted one.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        bool alpha = true,
        bool lower_star = false
    ) const;

    // Radii-free overload: plain (unweighted) Delaunay complex.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    // All multicolored simplices of the chosen complex. For full (weighted)
    // Delaunay triangulations every multicolored face is carried by a
    // multicolored tetrahedron, so only generating_tetrahedra is populated;
    // an alpha complex additionally yields free triangles and free edges.
    [[nodiscard]] MulticoloredSimplices collect_multicolored_simplices(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        bool alpha = true
    ) const;

    // Radii-free overload: plain (unweighted) Delaunay complex.
    [[nodiscard]] MulticoloredSimplices collect_multicolored_simplices(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

private:
    MulticoloredSimplices collect_multicolored_simplices_delaunay(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    MulticoloredSimplices collect_multicolored_simplices_weighted_delaunay(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    ) const;

    MulticoloredSimplices collect_multicolored_simplices_weighted_alpha(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    ) const;
};

[[nodiscard]] std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool alpha = true,
    bool lower_star = false
);

// Simplified 2-color surface: each vertex = midpoint of one cross-color atom pair.
// Requires exactly 2 distinct colors in the input.
[[nodiscard]] SimplifiedSurface compute_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool alpha = true
);

} // namespace delaunay_interfaces
