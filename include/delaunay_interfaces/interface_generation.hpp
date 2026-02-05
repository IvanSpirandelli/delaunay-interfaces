#pragma once

#include "types.hpp"
#include <memory>

namespace delaunay_interfaces {

class InterfaceGenerator {
public:
    InterfaceGenerator() = default;
    ~InterfaceGenerator() = default;

    // Main entry point
    InterfaceSurface compute_interface_surface(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii = {},
        bool weighted = true,
        bool alpha = true
    );

    // Get multicolored tetrahedra
    Tetrahedra get_multicolored_tetrahedra(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii = {},
        bool weighted = true,
        bool alpha = true
    );

private:
    // Delaunay/Alpha complex computation
    Tetrahedra get_multicolored_tetrahedra_delaunay(
        const Points& points,
        const ColorLabels& color_labels
    );

    Tetrahedra get_multicolored_tetrahedra_weighted_delaunay(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    );

    Tetrahedra get_multicolored_tetrahedra_weighted_alpha(
        const Points& points,
        const ColorLabels& color_labels,
        const Radii& radii
    );

    // Helper to check if tetrahedron is multicolored
    bool is_multicolored(const Tetrahedron& tet, const ColorLabels& color_labels) const;
};

// Barycentric subdivision functions
std::pair<Points, Filtration> get_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii = {},
    bool weighted = true,
    bool alpha = true
);

} // namespace delaunay_interfaces
