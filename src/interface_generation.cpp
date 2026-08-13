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

using UaVb = CGAL::Alpha_shape_vertex_base_3<Kernel, DtVb>;
using UaCb = CGAL::Alpha_shape_cell_base_3<Kernel, CGAL::Delaunay_triangulation_cell_base_3<Kernel>>;
using UaTds = CGAL::Triangulation_data_structure_3<UaVb, UaCb>;
using UaDt = CGAL::Delaunay_triangulation_3<Kernel, UaTds>;
using UnweightedAlphaShape = CGAL::Alpha_shape_3<UaDt>;
using IndexedPoint = std::pair<Kernel::Point_3, int>;

template <class VertexIndices>
static bool is_multicolored(const VertexIndices& indices, const ColorLabels& color_labels) {
    if (indices.size() < 2) return false;
    int first_color = color_labels[indices[0]];
    for (size_t i = 1; i < indices.size(); ++i) {
        if (color_labels[indices[i]] != first_color) return true;
    }
    return false;
}

// Tetrahedra of the complex whose vertices span >= 2 colors. keep_cell says
// which finite cells belong to the complex: all of them for a (weighted)
// Delaunay triangulation; for an alpha shape only cells not classified
// EXTERIOR, CGAL's label for simplices of the underlying triangulation that
// are not part of the alpha complex.
template <class Triangulation, class KeepCell>
static Tetrahedra collect_multicolored_tetrahedra(
    const Triangulation& t,
    const ColorLabels& color_labels,
    KeepCell keep_cell
) {
    Tetrahedra result;
    for (auto cit = t.finite_cells_begin(); cit != t.finite_cells_end(); ++cit) {
        if (!keep_cell(cit)) continue;

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

template <class Triangulation>
static Tetrahedra collect_multicolored_tetrahedra(
    const Triangulation& t,
    const ColorLabels& color_labels
) {
    return collect_multicolored_tetrahedra(t, color_labels, [](const auto&) { return true; });
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

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_delaunay(
    const Points& points,
    const ColorLabels& color_labels
) const {
    DelaunayTriangulation dt;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        auto vh = dt.insert(Kernel::Point_3(p.x(), p.y(), p.z()));
        vh->info() = static_cast<int>(i);
    }

    MulticoloredSimplices result;
    result.generating_tetrahedra = collect_multicolored_tetrahedra(dt, color_labels);
    return result;
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_weighted_delaunay(
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

    MulticoloredSimplices result;
    result.generating_tetrahedra = collect_multicolored_tetrahedra(rt, color_labels);
    return result;
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool alpha
) const {
    if (radii.empty()) {
        if (alpha) {
            throw std::invalid_argument(
                "alpha requires radii: pass per-point radii for the weighted "
                "alpha complex or a single radius for the uniform one");
        }
        return collect_multicolored_simplices_delaunay(points, color_labels);
    }
    return alpha
        ? collect_multicolored_simplices_weighted_alpha(points, color_labels, radii)
        : collect_multicolored_simplices_weighted_delaunay(points, color_labels, radii);
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return collect_multicolored_simplices(points, color_labels, Radii{}, false);
}

// Multicolored simplices of an alpha complex (weighted or unweighted).
// Free simplices are exactly the SINGULAR ones: in the alpha complex but
// not a face of any higher-dimensional simplex of the complex. (A
// multicolored face of an in-complex simplex is always covered by a
// collected multicolored simplex, since any simplex containing a
// multicolored face is itself multicolored.)
template <class AlphaShape>
static MulticoloredSimplices collect_from_alpha_shape(
    const AlphaShape& as,
    const ColorLabels& color_labels
) {
    MulticoloredSimplices result;
    result.generating_tetrahedra = collect_multicolored_tetrahedra(as, color_labels,
        [&as](const auto& cit) { return as.classify(cit) != AlphaShape::EXTERIOR; });

    // Free multicolored triangles.
    for (auto fit = as.finite_facets_begin(); fit != as.finite_facets_end(); ++fit) {
        if (as.classify(*fit) != AlphaShape::SINGULAR) continue;

        auto cell = fit->first;
        int face_idx = fit->second;
        FreeTriangle tri;
        int k = 0;
        for (int i = 0; i < 4; ++i) {
            if (i == face_idx) continue;
            tri[k++] = cell->vertex(i)->info();
        }

        if (is_multicolored(tri, color_labels)) {
            result.generating_free_triangles.push_back(tri);
        }
    }

    // Free multicolored edges.
    for (auto eit = as.finite_edges_begin(); eit != as.finite_edges_end(); ++eit) {
        if (as.classify(*eit) != AlphaShape::SINGULAR) continue;

        auto cell = eit->first;
        FreeEdge edge = {cell->vertex(eit->second)->info(), cell->vertex(eit->third)->info()};

        if (is_multicolored(edge, color_labels)) {
            result.generating_free_edges.push_back(edge);
        }
    }

    return result;
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) const {
    auto wpoints = make_indexed_weighted_points(points, radii);
    WeightedAlphaShape as(wpoints.begin(), wpoints.end(), 0, WeightedAlphaShape::GENERAL);
    return collect_from_alpha_shape(as, color_labels);
}

// Uniform-radius alpha complex: with equal weights the regular triangulation
// is the Delaunay triangulation and the alpha complex at alpha = 0 equals the
// unweighted one at alpha = radius^2, so the unweighted CGAL machinery
// (cheaper predicates, no hidden points) computes the same complex.
MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_uniform_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    double radius
) const {
    std::vector<IndexedPoint> ipoints;
    ipoints.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        ipoints.emplace_back(Kernel::Point_3(p.x(), p.y(), p.z()), static_cast<int>(i));
    }
    UnweightedAlphaShape as(ipoints.begin(), ipoints.end(), radius * radius,
                            UnweightedAlphaShape::GENERAL);
    return collect_from_alpha_shape(as, color_labels);
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool alpha
) const {
    // Uniform weights leave the Delaunay triangulation unchanged, so a
    // radius without alpha can only be a mistake.
    if (!alpha) {
        throw std::invalid_argument(
            "a radius has no effect on the plain Delaunay complex: omit it, "
            "or pass alpha=true for the alpha complex");
    }
    return collect_multicolored_simplices_uniform_alpha(points, color_labels, radius);
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool alpha,
    bool lower_star
) const {
    auto [vertices, filtration, simplices, vertex_atom_indices] = compute_barycentric_subdivision_and_filtration(
        points, color_labels, radii, alpha, lower_star
    );

    return InterfaceSurface{std::move(vertices), std::move(filtration), std::move(simplices),
                            std::move(vertex_atom_indices), !radii.empty(), alpha, lower_star};
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool alpha,
    bool lower_star
) const {
    auto [vertices, filtration, simplices, vertex_atom_indices] = compute_barycentric_subdivision_and_filtration(
        points, color_labels, radius, alpha, lower_star
    );

    return InterfaceSurface{std::move(vertices), std::move(filtration), std::move(simplices),
                            std::move(vertex_atom_indices), false, alpha, lower_star};
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return compute_interface_surface(points, color_labels, Radii{}, false, false);
}

} // namespace delaunay_interfaces
