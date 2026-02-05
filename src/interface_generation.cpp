#include "delaunay_interfaces/interface_generation.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <set>
#include <algorithm>

namespace delaunay_interfaces {

// CGAL type definitions - Delaunay triangulation
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using DelaunayVertexBase = CGAL::Triangulation_vertex_base_3<Kernel>;
using DelaunayCellBase = CGAL::Triangulation_cell_base_3<Kernel>;
using DelaunayDataStructure = CGAL::Triangulation_data_structure_3<DelaunayVertexBase, DelaunayCellBase>;
using DelaunayTriangulation = CGAL::Delaunay_triangulation_3<Kernel, DelaunayDataStructure>;

// Weighted Delaunay (Regular triangulation)
using WeightedPoint = Kernel::Weighted_point_3;
using RegularTriangulation = CGAL::Regular_triangulation_3<Kernel>;

// Alpha shapes
using AlphaVertexBase = CGAL::Triangulation_vertex_base_3<Kernel>;
using AlphaCellBase = CGAL::Alpha_shape_cell_base_3<Kernel, AlphaVertexBase>;
using AlphaDataStructure = CGAL::Triangulation_data_structure_3<AlphaVertexBase, AlphaCellBase>;
using AlphaTriangulation = CGAL::Delaunay_triangulation_3<Kernel, AlphaDataStructure>;
using AlphaShape = CGAL::Alpha_shape_3<AlphaTriangulation>;

bool InterfaceGenerator::is_multicolored(
    const Tetrahedron& tet,
    const ColorLabels& color_labels
) const {
    std::set<int> colors;
    for (int v : tet) {
        colors.insert(color_labels[v]);
    }
    return colors.size() >= 2;
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra_delaunay(
    const Points& points,
    const ColorLabels& color_labels
) {
    DelaunayTriangulation dt;
    std::map<DelaunayTriangulation::Vertex_handle, int> vertex_to_index;

    // Insert points and track indices
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        auto vh = dt.insert(Kernel::Point_3(p.x(), p.y(), p.z()));
        vertex_to_index[vh] = i;
    }

    Tetrahedra result;

    // Extract tetrahedra
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Tetrahedron tet;
        for (int i = 0; i < 4; ++i) {
            tet[i] = vertex_to_index[cit->vertex(i)];
        }

        if (is_multicolored(tet, color_labels)) {
            result.push_back(tet);
        }
    }

    return result;
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra_weighted_delaunay(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) {
    RegularTriangulation rt;
    std::map<RegularTriangulation::Vertex_handle, int> vertex_to_index;

    // Insert weighted points
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        double weight = radii[i] * radii[i]; // Weight is radius squared
        Kernel::Point_3 cgal_point(p.x(), p.y(), p.z());
        WeightedPoint wp(cgal_point, weight);

        auto vh = rt.insert(wp);
        if (vh != RegularTriangulation::Vertex_handle()) {
            vertex_to_index[vh] = i;
        }
    }

    Tetrahedra result;

    // Extract tetrahedra
    for (auto cit = rt.finite_cells_begin(); cit != rt.finite_cells_end(); ++cit) {
        if (rt.is_infinite(cit)) continue;

        Tetrahedron tet;
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            auto vh = cit->vertex(i);
            if (vertex_to_index.find(vh) == vertex_to_index.end()) {
                valid = false;
                break;
            }
            tet[i] = vertex_to_index[vh];
        }

        if (valid && is_multicolored(tet, color_labels)) {
            result.push_back(tet);
        }
    }

    return result;
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) {
    // For alpha shapes, we use regular triangulation and filter by alpha value
    RegularTriangulation rt;
    std::map<RegularTriangulation::Vertex_handle, int> vertex_to_index;

    // Insert weighted points
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        double weight = radii[i] * radii[i];
        Kernel::Point_3 cgal_point(p.x(), p.y(), p.z());
        WeightedPoint wp(cgal_point, weight);

        auto vh = rt.insert(wp);
        if (vh != RegularTriangulation::Vertex_handle()) {
            vertex_to_index[vh] = i;
        }
    }

    Tetrahedra result;

    // For weighted alpha complex with alpha=0, we include cells where the
    // weighted circumsphere has non-positive squared radius.
    // This is equivalent to the spheres (with given radii) already covering the cell.
    for (auto cit = rt.finite_cells_begin(); cit != rt.finite_cells_end(); ++cit) {
        if (rt.is_infinite(cit)) continue;

        Tetrahedron tet;
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            auto vh = cit->vertex(i);
            if (vertex_to_index.find(vh) == vertex_to_index.end()) {
                valid = false;
                break;
            }
            tet[i] = vertex_to_index[vh];
        }

        if (!valid || !is_multicolored(tet, color_labels)) {
            continue;
        }

        // Compute the critical alpha value for this cell
        // The power center (dual) of the cell and its power distance
        // A cell is in the alpha complex (alpha=0) if the squared distance
        // from the power center to any vertex (minus weight) is <= 0

        // Get the orthogonal center (power center) of the cell
        Kernel::Point_3 ortho_center = rt.dual(cit);

        // Check the squared power distance from ortho_center to any vertex
        // Power distance = squared_distance - weight
        auto v0 = cit->vertex(0);
        Kernel::Point_3 p0 = rt.point(v0).point();
        double w0 = rt.point(v0).weight();

        double sq_dist = CGAL::squared_distance(ortho_center, p0);
        double critical_value = sq_dist - w0;

        // Include cell if critical value <= 0 (matching Julia's fs <= 0)
        if (critical_value <= 0.0) {
            result.push_back(tet);
        }
    }

    return result;
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) {
    if (weighted) {
        if (alpha) {
            return get_multicolored_tetrahedra_weighted_alpha(points, color_labels, radii);
        } else {
            return get_multicolored_tetrahedra_weighted_delaunay(points, color_labels, radii);
        }
    } else {
        return get_multicolored_tetrahedra_delaunay(points, color_labels);
    }
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) {
    auto [vertices, filtration] = get_barycentric_subdivision_and_filtration(
        points, color_labels, radii, weighted, alpha
    );

    return InterfaceSurface{vertices, filtration, weighted, alpha};
}

} // namespace delaunay_interfaces
