#include "delaunay_interfaces/interface_generation.hpp"
#include "subdivision_driver.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_utils_3.h>
#include <algorithm>
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

// Alpha-path triangulations: the plain regular / Delaunay triangulation with
// the Fixed_alpha cell base, whose per-cell and per-facet classification
// slots our own classification pass fills in. CGAL::Fixed_alpha_shape_3 is no
// longer constructed — its initialize_alpha classifies every simplex up front
// (Gabriel predicate on nearly all facets before the cheap radius test, edge
// statuses for all edges in a global std::map) while the collection below
// needs only three queries. All classification semantics are replicated from
// CGAL 6.1 <CGAL/Fixed_alpha_shape_3.h>; line references below are into that
// header.
using AsCb = CGAL::Fixed_alpha_shape_cell_base_3<Kernel, CGAL::Regular_triangulation_cell_base_3<Kernel>>;
using AsTds = CGAL::Triangulation_data_structure_3<RtVb, AsCb>;
using AsRt = CGAL::Regular_triangulation_3<Kernel, AsTds>;
using IndexedWeightedPoint = std::pair<WeightedPoint, int>;

using UaCb = CGAL::Fixed_alpha_shape_cell_base_3<Kernel, CGAL::Delaunay_triangulation_cell_base_3<Kernel>>;
using UaTds = CGAL::Triangulation_data_structure_3<DtVb, UaCb>;
using UaDt = CGAL::Delaunay_triangulation_3<Kernel, UaTds>;
using IndexedPoint = std::pair<Kernel::Point_3, int>;

using CGAL::internal::EXTERIOR;
using CGAL::internal::INTERIOR;
using CGAL::internal::SINGULAR;

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

static std::vector<IndexedPoint> make_indexed_points(const Points& points) {
    std::vector<IndexedPoint> ipoints;
    ipoints.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        ipoints.emplace_back(Kernel::Point_3(p.x(), p.y(), p.z()), static_cast<int>(i));
    }
    return ipoints;
}

// The classification predicate, selected as in
// CGAL::internal::Simplex_classif_predicate (Fixed_alpha_shape_3.h:44-63):
// squared-radius comparison for the unweighted case, the weighted
// (power-distance) comparison for the regular case. Same kernel functors, so
// the filtered exact predicates give bit-identical answers.
static Kernel::Compare_squared_radius_3 alpha_predicate(const UaDt& t) {
    return t.geom_traits().compare_squared_radius_3_object();
}

static Kernel::Compare_weighted_squared_radius_3 alpha_predicate(const AsRt& t) {
    return t.geom_traits().compare_weighted_squared_radius_3_object();
}

// Classify every cell: EXTERIOR for infinite cells, otherwise INTERIOR iff
// the predicate on its four vertices vs alpha is != POSITIVE. Mirrors
// set_cell_status (Fixed_alpha_shape_3.h:816-820) calling
// is_gabriel_simplex_in_alpha_complex(Cell_handle) (lines 672-682): same
// functor, same vertex(0..3) argument order.
template <class Triangulation>
static void classify_cells(Triangulation& t, double alpha) {
    const auto in_complex = alpha_predicate(t);
    for (auto cit = t.all_cells_begin(); cit != t.all_cells_end(); ++cit) {
        if (t.is_infinite(cit)) {
            cit->set_classification_type(EXTERIOR);
        } else {
            cit->set_classification_type(
                in_complex(cit->vertex(0)->point(), cit->vertex(1)->point(),
                           cit->vertex(2)->point(), cit->vertex(3)->point(),
                           alpha) != CGAL::POSITIVE
                    ? INTERIOR : EXTERIOR);
        }
    }
}

// Multicolored simplices of an alpha complex (weighted or unweighted).
// Free simplices are exactly the SINGULAR ones: in the alpha complex but
// not a face of any higher-dimensional simplex of the complex. (A
// multicolored face of an in-complex simplex is always covered by a
// collected multicolored simplex, since any simplex containing a
// multicolored face is itself multicolored.)
//
// SINGULAR-ness is computed here instead of by Fixed_alpha_shape_3, with the
// identical exact predicates evaluated cheap-first — each status is a pure
// conjunction, so the order of evaluation cannot change the outcome:
//   facet SINGULAR (set_facet_classification_type, Fixed_alpha_shape_3.h:
//   835-873)  iff  no incident cell INTERIOR, is_Gabriel(f), and the radius
//   predicate on its vertex triple vs alpha is != POSITIVE;
//   edge SINGULAR (compute_edge_status, lines 887-916)  iff  every finite
//   incident facet EXTERIOR, is_Gabriel(c,i,j), and the radius predicate on
//   its two vertices vs alpha is != POSITIVE.
// "Every finite incident facet EXTERIOR" is tested as "no incident cell
// INTERIOR and no finite incident facet SINGULAR": a facet of an INTERIOR
// cell is INTERIOR or REGULAR, and once every incident cell is non-INTERIOR
// each finite incident facet is either SINGULAR or EXTERIOR.
template <class Triangulation>
static MulticoloredSimplices collect_from_alpha_complex(
    Triangulation& t,
    const ColorLabels& color_labels,
    double alpha
) {
    MulticoloredSimplices result;
    // Fixed_alpha_shape_3 only classifies 3-dimensional triangulations
    // (constructor guard, Fixed_alpha_shape_3.h:562-563).
    if (t.dimension() != 3) return result;

    classify_cells(t, alpha);
    const auto in_complex = alpha_predicate(t);

    result.tetrahedra = collect_multicolored_tetrahedra(t, color_labels,
        [](const auto& cit) { return cit->get_classification_type() != EXTERIOR; });

    // Free multicolored triangles. The cheap multicolor test runs first, then
    // the cell statuses, then the radius test; the Gabriel predicate only
    // when everything else says SINGULAR. The is-SINGULAR bit is cached in
    // both mirror facet slots (as set_facet_classification_type does) for the
    // edge pass below.
    for (auto fit = t.finite_facets_begin(); fit != t.finite_facets_end(); ++fit) {
        auto cell = fit->first;
        int face_idx = fit->second;
        FreeTriangle tri;
        int k = 0;
        for (int i = 0; i < 4; ++i) {
            if (i == face_idx) continue;
            tri[k++] = cell->vertex(i)->info();
        }
        // Canonical vertex order: which of the two incident cells reports the
        // facet is decided by raw Cell_handle pointer comparison (see the
        // note above the sorts at the end of this function), and that choice
        // also fixes the vertex_triple order.
        std::sort(tri.begin(), tri.end());
        if (!is_multicolored(tri, color_labels)) continue;

        // Radius arguments in vertex_triple_index order, as in
        // is_gabriel_simplex_in_alpha_complex(Cell_handle, int)
        // (Fixed_alpha_shape_3.h:684-693).
        auto neighbor = cell->neighbor(face_idx);
        const bool singular =
            cell->get_classification_type() != INTERIOR &&
            neighbor->get_classification_type() != INTERIOR &&
            in_complex(
                cell->vertex(CGAL::Triangulation_utils_3::vertex_triple_index(face_idx, 0))->point(),
                cell->vertex(CGAL::Triangulation_utils_3::vertex_triple_index(face_idx, 1))->point(),
                cell->vertex(CGAL::Triangulation_utils_3::vertex_triple_index(face_idx, 2))->point(),
                alpha) != CGAL::POSITIVE &&
            t.is_Gabriel(*fit);

        const auto status = singular ? SINGULAR : EXTERIOR;
        cell->set_facet_classification_type(face_idx, status);
        neighbor->set_facet_classification_type(neighbor->index(cell), status);

        if (singular) result.free_triangles.push_back(tri);
    }

    // Free multicolored edges, multicolor test first as above. Every facet
    // incident to a multicolored edge contains that edge's two differently
    // colored vertices, hence is multicolored, hence its is-SINGULAR bit was
    // cached by the facet pass.
    for (auto eit = t.finite_edges_begin(); eit != t.finite_edges_end(); ++eit) {
        auto cell = eit->first;
        FreeEdge edge = {cell->vertex(eit->second)->info(), cell->vertex(eit->third)->info()};
        std::sort(edge.begin(), edge.end()); // canonical, as for the facets above
        if (!is_multicolored(edge, color_labels)) continue;

        auto ccirc = t.incident_cells(cell, eit->second, eit->third);
        auto cdone = ccirc;
        bool candidate = true;
        do {
            if (ccirc->get_classification_type() == INTERIOR) {
                candidate = false;
                break;
            }
        } while (++ccirc != cdone);
        if (!candidate) continue;

        // Infinite incident facets are skipped exactly as in
        // compute_edge_status (Fixed_alpha_shape_3.h:899).
        auto fcirc = t.incident_facets(cell, eit->second, eit->third);
        auto fdone = fcirc;
        do {
            if (!t.is_infinite(*fcirc) &&
                fcirc->first->get_facet_classification_type(fcirc->second) == SINGULAR) {
                candidate = false;
                break;
            }
        } while (++fcirc != fdone);
        if (!candidate) continue;

        // Radius arguments as in is_gabriel_simplex_in_alpha_complex
        // (Cell_handle, int, int) (Fixed_alpha_shape_3.h:700-708); Gabriel
        // call as in line 912.
        if (in_complex(cell->vertex(eit->second)->point(),
                       cell->vertex(eit->third)->point(),
                       alpha) != CGAL::POSITIVE &&
            t.is_Gabriel(cell, eit->second, eit->third)) {
            result.free_edges.push_back(edge);
        }
    }

    // The two free lists are emitted in triangulation traversal order, which
    // is not reproducible: CGAL's finite facet / edge iterators report each
    // facet (resp. edge) from the incident cell with the *smaller pointer*
    // (Triangulation_ds_iterators_3.h:83 and :232). Cell_handle comparison
    // falls back to No_time_stamp::less — raw addresses — because
    // Fixed_alpha_shape_cell_base_3 declares no Has_timestamp tag, so the
    // Compact_container block layout, and with it the emission order, varies
    // between two identical calls in one process. Sorting here makes the
    // lists a function of the atom indices alone. The tetrahedra are left
    // untouched: finite_cells_begin() is plain container order, with no
    // pointer comparison, so it is already reproducible — and since the
    // subdivision feeds tetrahedra first, every vertex id they create keeps
    // its current value.
    std::sort(result.free_triangles.begin(), result.free_triangles.end());
    std::sort(result.free_edges.begin(), result.free_edges.end());
    return result;
}

// Both constructions use CGAL's range insert on (point, index) pairs, which
// Hilbert-sorts the input first — much faster than inserting unordered
// points one at a time, and it attaches the info (and drops hidden points'
// info, for the regular triangulation) itself.
MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_delaunay(
    const Points& points,
    const ColorLabels& color_labels
) const {
    auto ipoints = make_indexed_points(points);
    DelaunayTriangulation dt(ipoints.begin(), ipoints.end());

    MulticoloredSimplices result;
    result.tetrahedra = collect_multicolored_tetrahedra(dt, color_labels);
    return result;
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_weighted_delaunay(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) const {
    auto wpoints = make_indexed_weighted_points(points, radii);
    RegularTriangulation rt(wpoints.begin(), wpoints.end());

    MulticoloredSimplices result;
    result.tetrahedra = collect_multicolored_tetrahedra(rt, color_labels);
    return result;
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) const {
    auto wpoints = make_indexed_weighted_points(points, radii);
    AsRt rt(wpoints.begin(), wpoints.end());
    return collect_from_alpha_complex(rt, color_labels, 0.0);
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
    auto ipoints = make_indexed_points(points);
    UaDt dt(ipoints.begin(), ipoints.end());
    return collect_from_alpha_complex(dt, color_labels, radius * radius);
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha
) const {
    const bool use_alpha = alpha.value_or(!radii.empty());
    if (radii.empty()) {
        if (use_alpha) {
            throw std::invalid_argument(
                "alpha requires radii: pass per-point radii for the weighted "
                "alpha complex or a single radius for the uniform one");
        }
        return collect_multicolored_simplices_delaunay(points, color_labels);
    }
    return use_alpha
        ? collect_multicolored_simplices_weighted_alpha(points, color_labels, radii)
        : collect_multicolored_simplices_weighted_delaunay(points, color_labels, radii);
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels,
    double radius
) const {
    return collect_multicolored_simplices_uniform_alpha(points, color_labels, radius);
}

MulticoloredSimplices InterfaceGenerator::collect_multicolored_simplices(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return collect_multicolored_simplices(points, color_labels, Radii{}, false);
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha,
    bool lower_star
) const {
    const bool use_alpha = alpha.value_or(!radii.empty());
    // Flat pipeline: the surface keeps the compact filtration form and only
    // materializes the vector-of-vectors view on demand.
    auto [vertices, filtration, simplices, vertex_atom_indices] = detail::compute_subdivision_flat(
        points, color_labels, radii, use_alpha, lower_star
    );

    return InterfaceSurface{std::move(vertices), std::move(filtration), std::move(simplices),
                            std::move(vertex_atom_indices), use_alpha, lower_star};
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool lower_star
) const {
    auto [vertices, filtration, simplices, vertex_atom_indices] = detail::compute_subdivision_flat(
        points, color_labels, radius, lower_star
    );

    return InterfaceSurface{std::move(vertices), std::move(filtration), std::move(simplices),
                            std::move(vertex_atom_indices), true, lower_star};
}

InterfaceSurface InterfaceGenerator::compute_interface_surface(
    const Points& points,
    const ColorLabels& color_labels
) const {
    return compute_interface_surface(points, color_labels, Radii{}, false, false);
}

} // namespace delaunay_interfaces
