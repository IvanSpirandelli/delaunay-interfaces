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

// Helper: check if a set of vertex indices is multicolored
static bool is_multicolored_set(const std::vector<int>& indices, const ColorLabels& color_labels) {
    if (indices.size() < 2) return false;
    int first_color = color_labels[indices[0]];
    for (size_t i = 1; i < indices.size(); ++i) {
        if (color_labels[indices[i]] != first_color) return true;
    }
    return false;
}

// Helper: make a sorted vertex index vector from a set of indices
static std::vector<int> make_sorted_key(std::initializer_list<int> indices) {
    std::vector<int> key(indices);
    std::sort(key.begin(), key.end());
    return key;
}

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

// Helper: compute alpha value (critical squared power radius) for a cell
static double compute_cell_alpha(
    const RegularTriangulation& rt,
    RegularTriangulation::Cell_handle cell
) {
    Kernel::Point_3 ortho_center = rt.dual(cell);
    auto v0 = cell->vertex(0);
    Kernel::Point_3 p0 = rt.point(v0).point();
    double w0 = rt.point(v0).weight();
    double sq_dist = CGAL::squared_distance(ortho_center, p0);
    return sq_dist - w0;
}

// Helper: compute alpha value for a facet (triangle) using weighted circumcenter
static double compute_facet_alpha(
    const WeightedPoint& wp0,
    const WeightedPoint& wp1,
    const WeightedPoint& wp2
) {
    // The weighted circumcenter is the point in the plane of the 3 points
    // equidistant (in power distance) from all three weighted points.
    Kernel::Construct_weighted_circumcenter_3 wcc;
    Kernel::Point_3 center = wcc(wp0, wp1, wp2);

    double sq_dist = CGAL::squared_distance(center, wp0.point());
    return sq_dist - wp0.weight();
}

// Helper: check if a point lies inside a triangle (in 3D, assumes coplanar)
static bool point_in_triangle(
    const Kernel::Point_3& p,
    const Kernel::Point_3& a,
    const Kernel::Point_3& b,
    const Kernel::Point_3& c
) {
    // Use barycentric coordinates
    auto v0 = b - a;
    auto v1 = c - a;
    auto v2 = p - a;

    double d00 = v0 * v0;
    double d01 = v0 * v1;
    double d11 = v1 * v1;
    double d20 = v2 * v0;
    double d21 = v2 * v1;

    double denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1e-15) return false; // degenerate triangle

    double u = (d11 * d20 - d01 * d21) / denom;
    double v = (d00 * d21 - d01 * d20) / denom;

    return (u >= -1e-10) && (v >= -1e-10) && (u + v <= 1.0 + 1e-10);
}

// Helper: compute alpha value for an edge using weighted circumcenter
static double compute_edge_alpha(
    const WeightedPoint& wp0,
    const WeightedPoint& wp1
) {
    Kernel::Construct_weighted_circumcenter_3 wcc;
    Kernel::Point_3 center = wcc(wp0, wp1);

    double sq_dist = CGAL::squared_distance(center, wp0.point());
    return sq_dist - wp0.weight();
}

// Helper: check if weighted circumcenter lies on the edge segment
static bool edge_is_gabriel(
    const WeightedPoint& wp0,
    const WeightedPoint& wp1
) {
    Kernel::Construct_weighted_circumcenter_3 wcc;
    Kernel::Point_3 center = wcc(wp0, wp1);

    // Check if center lies between p0 and p1 using parameter t
    auto p0 = wp0.point();
    auto p1 = wp1.point();
    auto edge_vec = p1 - p0;
    auto to_center = center - p0;

    double edge_sq = edge_vec * edge_vec;
    if (edge_sq < 1e-15) return false; // degenerate edge

    double t = (to_center * edge_vec) / edge_sq;
    return (t >= -1e-10) && (t <= 1.0 + 1e-10);
}

MulticoloredSimplices InterfaceGenerator::get_multicolored_simplices_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) {
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

    MulticoloredSimplices result;

    // --- Pass 1: Extract alpha cells and multicolored alpha tetrahedra ---
    std::set<RegularTriangulation::Cell_handle> alpha_cells;

    for (auto cit = rt.finite_cells_begin(); cit != rt.finite_cells_end(); ++cit) {
        if (rt.is_infinite(cit)) continue;

        // Check all vertices are mapped
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            if (vertex_to_index.find(cit->vertex(i)) == vertex_to_index.end()) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        double alpha_val = compute_cell_alpha(rt, cit);
        if (alpha_val <= 0.0) {
            alpha_cells.insert(cit);

            // Build tetrahedron
            Tetrahedron tet;
            for (int i = 0; i < 4; ++i) {
                tet[i] = vertex_to_index[cit->vertex(i)];
            }

            if (is_multicolored(tet, color_labels)) {
                result.tetrahedra.push_back(tet);
            }
        }
    }

    // Build set of triangle faces covered by multicolored alpha tetrahedra
    std::set<std::vector<int>> covered_tris;
    // Build set of edges covered by multicolored alpha tetrahedra
    std::set<std::vector<int>> covered_edges;

    for (const auto& tet : result.tetrahedra) {
        // 4 triangle faces
        for (int skip = 0; skip < 4; ++skip) {
            std::vector<int> tri;
            for (int i = 0; i < 4; ++i) {
                if (i != skip) tri.push_back(tet[i]);
            }
            std::sort(tri.begin(), tri.end());
            covered_tris.insert(tri);
        }
        // 6 edges
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                covered_edges.insert(make_sorted_key({tet[i], tet[j]}));
            }
        }
    }

    // --- Pass 2: Extract free multicolored facets ---
    for (auto fit = rt.finite_facets_begin(); fit != rt.finite_facets_end(); ++fit) {
        auto cell = fit->first;
        int face_idx = fit->second;

        // Get the 3 vertices of this facet
        std::vector<int> tri_indices;
        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            if (i == face_idx) continue;
            auto vh = cell->vertex(i);
            auto it = vertex_to_index.find(vh);
            if (it == vertex_to_index.end()) { valid = false; break; }
            tri_indices.push_back(it->second);
        }
        if (!valid || tri_indices.size() != 3) continue;

        // Check if multicolored
        if (!is_multicolored_set(tri_indices, color_labels)) continue;

        // Check if already covered by a multicolored alpha tetrahedron
        std::vector<int> sorted_tri = tri_indices;
        std::sort(sorted_tri.begin(), sorted_tri.end());
        if (covered_tris.count(sorted_tri)) continue;

        // Check if in the alpha complex:
        // (a) Adjacent cell is in alpha complex
        bool in_alpha = alpha_cells.count(cell) > 0;
        if (!in_alpha) {
            auto neighbor = cell->neighbor(face_idx);
            if (!rt.is_infinite(neighbor)) {
                in_alpha = alpha_cells.count(neighbor) > 0;
            }
        }

        // (b) Gabriel facet with alpha <= 0
        if (!in_alpha) {
            auto vh0 = cell->vertex((face_idx + 1) % 4);
            auto vh1 = cell->vertex((face_idx + 2) % 4);
            auto vh2 = cell->vertex((face_idx + 3) % 4);

            WeightedPoint wp0 = rt.point(vh0);
            WeightedPoint wp1 = rt.point(vh1);
            WeightedPoint wp2 = rt.point(vh2);

            double facet_alpha = compute_facet_alpha(wp0, wp1, wp2);
            if (facet_alpha <= 0.0) {
                // Gabriel check: circumcenter inside triangle
                Kernel::Construct_weighted_circumcenter_3 wcc;
                Kernel::Point_3 center = wcc(wp0, wp1, wp2);
                if (point_in_triangle(center, wp0.point(), wp1.point(), wp2.point())) {
                    in_alpha = true;
                }
            }
        }

        if (in_alpha) {
            result.free_triangles.push_back(tri_indices);
            // Add edges of this free triangle to covered_edges
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    covered_edges.insert(make_sorted_key({tri_indices[i], tri_indices[j]}));
                }
            }
        }
    }

    // --- Pass 3: Extract free multicolored edges ---
    for (auto eit = rt.finite_edges_begin(); eit != rt.finite_edges_end(); ++eit) {
        auto cell = eit->first;
        int i = eit->second;
        int j = eit->third;

        auto vh0 = cell->vertex(i);
        auto vh1 = cell->vertex(j);

        auto it0 = vertex_to_index.find(vh0);
        auto it1 = vertex_to_index.find(vh1);
        if (it0 == vertex_to_index.end() || it1 == vertex_to_index.end()) continue;

        int idx0 = it0->second;
        int idx1 = it1->second;

        // Check multicolored
        if (color_labels[idx0] == color_labels[idx1]) continue;

        // Check if already covered
        auto sorted_edge = make_sorted_key({idx0, idx1});
        if (covered_edges.count(sorted_edge)) continue;

        // Check if in the alpha complex:
        // (a) Any incident cell is in alpha complex
        bool in_alpha = false;
        auto circ = rt.incident_cells(*eit);
        auto start = circ;
        do {
            if (!rt.is_infinite(circ)) {
                if (alpha_cells.count(circ) > 0) {
                    in_alpha = true;
                    break;
                }
            }
            ++circ;
        } while (circ != start);

        // (b) Gabriel edge with alpha <= 0
        if (!in_alpha) {
            WeightedPoint wp0 = rt.point(vh0);
            WeightedPoint wp1 = rt.point(vh1);

            double edge_alpha = compute_edge_alpha(wp0, wp1);
            if (edge_alpha <= 0.0 && edge_is_gabriel(wp0, wp1)) {
                in_alpha = true;
            }
        }

        if (in_alpha) {
            result.free_edges.push_back({idx0, idx1});
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
) {
    auto [vertices, filtration, simplices] = get_barycentric_subdivision_and_filtration(
        points, color_labels, radii, weighted, alpha, lower_star
    );

    return InterfaceSurface{vertices, filtration, std::move(simplices), weighted, alpha, lower_star};
}

} // namespace delaunay_interfaces
