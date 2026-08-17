#pragma once

#include "types.hpp"
#include <optional>

namespace delaunay_interfaces {

class InterfaceGenerator {
public:
    // Weightedness is determined by the radii: empty radii build the plain
    // Delaunay complex, non-empty radii the weighted one. Unset alpha is
    // derived from the radii: alpha complex iff radii are given. Explicitly
    // passing alpha=true with empty radii throws (alpha requires radii);
    // alpha=false with radii selects the weighted Delaunay complex.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        std::optional<bool> alpha = std::nullopt,
        bool lower_star = false
    ) const;

    // Uniform-radius overload: builds the alpha complex with parameter
    // radius^2 via CGAL's unweighted alpha shape. Always an alpha complex:
    // a radius has no effect on the plain Delaunay complex, so there is no
    // alpha parameter here; use the radii-free overload for plain Delaunay.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels,
        double radius,
        bool lower_star = false
    ) const;

    // Radii-free overload: plain (unweighted) Delaunay complex.
    [[nodiscard]] InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels
    ) const;

    // All multicolored simplices of the chosen complex. For full (weighted)
    // Delaunay triangulations every multicolored face is carried by a
    // multicolored tetrahedron, so only tetrahedra is populated; an alpha
    // complex additionally yields free triangles and free edges.
    [[nodiscard]] MulticoloredSimplices collect_multicolored_simplices(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii,
        std::optional<bool> alpha = std::nullopt
    ) const;

    // Uniform-radius overload; see compute_interface_surface.
    [[nodiscard]] MulticoloredSimplices collect_multicolored_simplices(
        const Points& points,
        const ColorLabels& color_labels,
        double radius
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

    MulticoloredSimplices collect_multicolored_simplices_uniform_alpha(
        const Points& points,
        const ColorLabels& color_labels,
        double radius
    ) const;
};

[[nodiscard]] std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha = std::nullopt,
    bool lower_star = false
);

// Uniform-radius overload; see InterfaceGenerator::compute_interface_surface.
[[nodiscard]] std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool lower_star = false
);

// Simplified 2-color surface: each vertex = midpoint of one cross-color atom pair.
// Requires exactly 2 distinct colors in the input.
[[nodiscard]] SimplifiedSurface compute_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha = std::nullopt
);

// Uniform-radius overload; see InterfaceGenerator::compute_interface_surface.
[[nodiscard]] SimplifiedSurface compute_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    double radius
);

} // namespace delaunay_interfaces
