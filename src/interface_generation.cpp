#include "delaunay_interfaces/interface_generation.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <stdexcept>

namespace delaunay_interfaces {

// CGAL aliases. All triangulations carry the input point index as vertex info.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using DtVb = CGAL::Triangulation_vertex_base_with_info_3<int, Kernel>;
using DtTds = CGAL::Triangulation_data_structure_3<DtVb, CGAL::Delaunay_triangulation_cell_base_3<Kernel>>;
using DelaunayTriangulation = CGAL::Delaunay_triangulation_3<Kernel, DtTds>;

using WeightedPoint = Kernel::Weighted_point_3;
using RtVb = CGAL::Triangulation_vertex_base_with_info_3<int, Kernel, CGAL::Regular_triangulation_vertex_base_3<Kernel>>;
using RtTds = CGAL::Triangulation_data_structure_3<RtVb, CGAL::Regular_triangulation_cell_base_3<Kernel>>;
using RegularTriangulation = CGAL::Regular_triangulation_3<Kernel, RtTds>;

using AsVb = CGAL::Alpha_shape_vertex_base_3<Kernel, RtVb>;
using AsCb = CGAL::Alpha_shape_cell_base_3<Kernel, CGAL::Regular_triangulation_cell_base_3<Kernel>>;
using AsTds = CGAL::Triangulation_data_structure_3<AsVb, AsCb>;
using AsRt = CGAL::Regular_triangulation_3<Kernel, AsTds>;
using WeightedAlphaShape = CGAL::Alpha_shape_3<AsRt>;
using IndexedWeightedPoint = std::pair<WeightedPoint, int>;

static bool is_multicolored_set(const std::vector<int>& indices, const ColorLabels& color_labels) {
    if (indices.size() < 2) return false;
    int first_color = color_labels[indices[0]];
    for (size_t i = 1; i < indices.size(); ++i) {
        if (color_labels[indices[i]] != first_color) return true;
    }
    return false;
}

static bool is_multicolored_tet(const Tetrahedron& tet, const ColorLabels& color_labels) {
    int first_color = color_labels[tet[0]];
    return color_labels[tet[1]] != first_color
        || color_labels[tet[2]] != first_color
        || color_labels[tet[3]] != first_color;
}

static std::vector<IndexedWeightedPoint> make_indexed_weighted_points(
    const Points& points,
    const Radii& radii
) {
    std::vector<IndexedWeightedPoint> wpoints;
    wpoints.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        double weight = radii[i] * radii[i];
        wpoints.emplace_back(
            WeightedPoint(Kernel::Point_3(p.x(), p.y(), p.z()), weight),
            static_cast<int>(i)
        );
    }
    return wpoints;
}

// Tetrahedra of the alpha shape (non-EXTERIOR cells) whose vertices span >= 2 colors.
static Tetrahedra collect_multicolored_alpha_tetrahedra(
    const WeightedAlphaShape& as,
    const ColorLabels& color_labels
) {
    Tetrahedra result;
    for (auto cit = as.finite_cells_begin(); cit != as.finite_cells_end(); ++cit) {
        if (as.classify(cit) == WeightedAlphaShape::EXTERIOR) continue;

        Tetrahedron tet;
        for (int i = 0; i < 4; ++i) {
            tet[i] = cit->vertex(i)->info();
        }
        if (is_multicolored_tet(tet, color_labels)) {
            result.push_back(tet);
        }
    }
    return result;
}

bool InterfaceGenerator::is_multicolored(
    const Tetrahedron& tet,
    const ColorLabels& color_labels
) const {
    return is_multicolored_tet(tet, color_labels);
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra_delaunay(
    const Points& points,
    const ColorLabels& color_labels
) const {
    DelaunayTriangulation dt;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        auto vh = dt.insert(Kernel::Point_3(p.x(), p.y(), p.z()));
        vh->info() = static_cast<int>(i);
    }

    Tetrahedra result;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Tetrahedron tet;
        for (int i = 0; i < 4; ++i) {
            tet[i] = cit->vertex(i)->info();
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
) const {
    RegularTriangulation rt;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        WeightedPoint wp(Kernel::Point_3(p.x(), p.y(), p.z()), radii[i] * radii[i]);
        auto vh = rt.insert(wp);
        // Points hidden by the regular triangulation yield a null handle.
        if (vh != RegularTriangulation::Vertex_handle()) {
            vh->info() = static_cast<int>(i);
        }
    }

    Tetrahedra result;
    for (auto cit = rt.finite_cells_begin(); cit != rt.finite_cells_end(); ++cit) {
        Tetrahedron tet;
        for (int i = 0; i < 4; ++i) {
            tet[i] = cit->vertex(i)->info();
        }
        if (is_multicolored(tet, color_labels)) {
            result.push_back(tet);
        }
    }
    return result;
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) const {
    auto wpoints = make_indexed_weighted_points(points, radii);
    WeightedAlphaShape as(wpoints.begin(), wpoints.end(), 0, WeightedAlphaShape::GENERAL);
    return collect_multicolored_alpha_tetrahedra(as, color_labels);
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) const {
    if (!weighted && alpha) {
        throw std::invalid_argument("alpha complexes require weighted=true (unweighted alpha is not implemented)");
    }
    if (weighted) {
        return alpha
            ? get_multicolored_tetrahedra_weighted_alpha(points, color_labels, radii)
            : get_multicolored_tetrahedra_weighted_delaunay(points, color_labels, radii);
    }
    return get_multicolored_tetrahedra_delaunay(points, color_labels);
}

Tetrahedra InterfaceGenerator::get_multicolored_tetrahedra(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return get_multicolored_tetrahedra(points, color_labels, {}, false, false);
}

MulticoloredSimplices InterfaceGenerator::get_multicolored_simplices_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) const {
    auto wpoints = make_indexed_weighted_points(points, radii);
    WeightedAlphaShape as(wpoints.begin(), wpoints.end(), 0, WeightedAlphaShape::GENERAL);

    MulticoloredSimplices result;
    result.generating_tetrahedra = collect_multicolored_alpha_tetrahedra(as, color_labels);

    // Free simplices are exactly the SINGULAR ones: in the alpha complex but
    // not a face of any higher-dimensional simplex of the complex. (A
    // multicolored face of an in-complex simplex is always covered by a
    // collected multicolored simplex, since any simplex containing a
    // multicolored face is itself multicolored.)

    // Free multicolored triangles.
    for (auto fit = as.finite_facets_begin(); fit != as.finite_facets_end(); ++fit) {
        if (as.classify(*fit) != WeightedAlphaShape::SINGULAR) continue;

        auto cell = fit->first;
        int face_idx = fit->second;
        std::vector<int> tri_indices;
        for (int i = 0; i < 4; ++i) {
            if (i == face_idx) continue;
            tri_indices.push_back(cell->vertex(i)->info());
        }

        if (is_multicolored_set(tri_indices, color_labels)) {
            result.generating_free_triangles.push_back(tri_indices);
        }
    }

    // Free multicolored edges.
    for (auto eit = as.finite_edges_begin(); eit != as.finite_edges_end(); ++eit) {
        if (as.classify(*eit) != WeightedAlphaShape::SINGULAR) continue;

        auto cell = eit->first;
        int idx0 = cell->vertex(eit->second)->info();
        int idx1 = cell->vertex(eit->third)->info();

        if (color_labels[idx0] != color_labels[idx1]) {
            result.generating_free_edges.push_back({idx0, idx1});
        }
    }

    return result;
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha,
    bool lower_star
) const {
    auto [vertices, filtration, simplices, vertex_atom_indices] = get_barycentric_subdivision_and_filtration(
        points, color_labels, radii, weighted, alpha, lower_star
    );

    return InterfaceSurface{std::move(vertices), std::move(filtration), std::move(simplices),
                            std::move(vertex_atom_indices), weighted, alpha, lower_star};
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return compute_interface_surface(points, color_labels, {}, false, false, false);
}

} // namespace delaunay_interfaces
