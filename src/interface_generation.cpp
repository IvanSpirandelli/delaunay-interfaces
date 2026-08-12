#include "delaunay_interfaces/interface_generation.hpp"
#include "delaunay_interfaces/simplified_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
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

// Weighted alpha shapes (Regular_triangulation_3 based, with vertex info for index tracking)
using RtVb = CGAL::Regular_triangulation_vertex_base_3<Kernel>;
using RtVbInfo = CGAL::Triangulation_vertex_base_with_info_3<int, Kernel, RtVb>;
using AsVb = CGAL::Alpha_shape_vertex_base_3<Kernel, RtVbInfo>;
using RtCb = CGAL::Regular_triangulation_cell_base_3<Kernel>;
using AsCb = CGAL::Alpha_shape_cell_base_3<Kernel, RtCb>;
using AsTds = CGAL::Triangulation_data_structure_3<AsVb, AsCb>;
using AsRt = CGAL::Regular_triangulation_3<Kernel, AsTds>;
using WeightedAlphaShape = CGAL::Alpha_shape_3<AsRt>;
using IndexedWeightedPoint = std::pair<WeightedPoint, int>;

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

// Build indexed weighted points for alpha shape construction
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
    auto wpoints = make_indexed_weighted_points(points, radii);
    WeightedAlphaShape as(wpoints.begin(), wpoints.end(), 0, WeightedAlphaShape::GENERAL);

    Tetrahedra result;

    for (auto cit = as.finite_cells_begin(); cit != as.finite_cells_end(); ++cit) {
        if (as.classify(cit) == WeightedAlphaShape::EXTERIOR) continue;

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

MulticoloredSimplices InterfaceGenerator::get_multicolored_simplices_weighted_alpha(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii
) {
    auto wpoints = make_indexed_weighted_points(points, radii);
    WeightedAlphaShape as(wpoints.begin(), wpoints.end(), 0, WeightedAlphaShape::GENERAL);

    MulticoloredSimplices result;

    // --- Pass 1: Extract multicolored alpha tetrahedra ---
    for (auto cit = as.finite_cells_begin(); cit != as.finite_cells_end(); ++cit) {
        if (as.classify(cit) == WeightedAlphaShape::EXTERIOR) continue;

        Tetrahedron tet;
        for (int i = 0; i < 4; ++i) {
            tet[i] = cit->vertex(i)->info();
        }

        if (is_multicolored(tet, color_labels)) {
            result.generating_tetrahedra.push_back(tet);
        }
    }

    // Build set of triangle faces covered by multicolored alpha tetrahedra
    std::set<std::vector<int>> covered_tris;
    // Build set of edges covered by multicolored alpha tetrahedra or free triangles
    std::set<std::vector<int>> covered_edges;

    for (const auto& tet : result.generating_tetrahedra) {
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
    for (auto fit = as.finite_facets_begin(); fit != as.finite_facets_end(); ++fit) {
        if (as.classify(*fit) == WeightedAlphaShape::EXTERIOR) continue;

        auto cell = fit->first;
        int face_idx = fit->second;

        // Get the 3 vertex indices of this facet
        std::vector<int> tri_indices;
        for (int i = 0; i < 4; ++i) {
            if (i == face_idx) continue;
            tri_indices.push_back(cell->vertex(i)->info());
        }

        // Check if multicolored
        if (!is_multicolored_set(tri_indices, color_labels)) continue;

        // Check if already covered by a multicolored alpha tetrahedron
        std::vector<int> sorted_tri = tri_indices;
        std::sort(sorted_tri.begin(), sorted_tri.end());
        if (covered_tris.count(sorted_tri)) continue;

        result.generating_free_triangles.push_back(tri_indices);
        // Add edges of this free triangle to covered_edges
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                covered_edges.insert(make_sorted_key({tri_indices[i], tri_indices[j]}));
            }
        }
    }

    // --- Pass 3: Extract free multicolored edges ---
    for (auto eit = as.finite_edges_begin(); eit != as.finite_edges_end(); ++eit) {
        if (as.classify(*eit) == WeightedAlphaShape::EXTERIOR) continue;

        auto cell = eit->first;
        int idx0 = cell->vertex(eit->second)->info();
        int idx1 = cell->vertex(eit->third)->info();

        // Check multicolored
        if (color_labels[idx0] == color_labels[idx1]) continue;

        // Check if already covered
        auto sorted_edge = make_sorted_key({idx0, idx1});
        if (covered_edges.count(sorted_edge)) continue;

        result.generating_free_edges.push_back({idx0, idx1});
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
    auto [vertices, filtration, simplices, vertex_atom_indices] = get_barycentric_subdivision_and_filtration(
        points, color_labels, radii, weighted, alpha, lower_star
    );

    return InterfaceSurface{vertices, filtration, std::move(simplices), std::move(vertex_atom_indices), weighted, alpha, lower_star};
}

SimplifiedSurface get_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) {
    if (points.size() != color_labels.size()) {
        throw std::invalid_argument("Each point must have a corresponding color_label");
    }
    if (weighted && radii.size() != points.size()) {
        throw std::invalid_argument("Each point must have an assigned radius for weighted complexes");
    }

    InterfaceGenerator generator;
    SimplifiedSubdivision subdivision(points, color_labels);
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

        mc_simplices.generating_tetrahedra = tetrahedra;

        for (const auto& tet : tetrahedra) {
            subdivision.process_tetrahedron(tet);
        }
    }

    return SimplifiedSurface{
        subdivision.get_vertices(),
        subdivision.get_triangles(),
        subdivision.get_quads(),
        subdivision.get_edges(),
        subdivision.get_vertex_atom_indices(),
        subdivision.get_vertex_filtration(),
        std::move(mc_simplices),
        weighted,
        alpha
    };
}

} // namespace delaunay_interfaces
