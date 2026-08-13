#include <iostream>
#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <random>
#include <delaunay_interfaces/interface_generation.hpp>
#include <delaunay_interfaces/simplified_subdivision.hpp>
#include <delaunay_interfaces/barycentric_subdivision.hpp>
#include <delaunay_interfaces/chromatic_partitioning.hpp>

using namespace delaunay_interfaces;

void test_simplified_22_tetrahedron() {
    std::cout << "Test: Simplified 2-2 Tetrahedron\n";

    Points points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // color 1
        {0.5, 1.0, 0.0}, {0.5, 0.5, 1.0}    // color 2
    };
    ColorLabels colors = {1, 1, 2, 2};

    SimplifiedSubdivision sub(points, colors);
    sub.process_tetrahedron({0, 1, 2, 3});

    // 2-2: 4 cross-color edges → 4 vertices, 1 quad
    assert(sub.get_midpoints().size() == 4);
    assert(sub.get_quads().size() == 1);
    assert(sub.get_triangles().size() == 0);
    assert(sub.get_edges().size() == 0);

    // Every vertex must have exactly 2 atom indices
    auto vai = sub.get_vertex_atom_indices();
    assert(vai.size() == 4);
    for (const auto& atoms : vai) {
        assert(atoms.size() == 2);
        assert(colors[atoms[0]] != colors[atoms[1]]);
    }

    // Filtration values should be positive distances
    auto filt = sub.get_vertex_filtration();
    assert(filt.size() == 4);
    for (double v : filt) {
        assert(v > 0.0);
    }

    std::cout << "  4 vertices, 1 quad - PASS\n";
}

void test_simplified_31_tetrahedron() {
    std::cout << "Test: Simplified 3-1 Tetrahedron\n";

    Points points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 1.0, 0.0},  // color 1
        {0.5, 0.5, 1.0}                                         // color 2
    };
    ColorLabels colors = {1, 1, 1, 2};

    SimplifiedSubdivision sub(points, colors);
    sub.process_tetrahedron({0, 1, 2, 3});

    // 3-1: 3 cross-color edges → 3 vertices, 1 triangle
    assert(sub.get_midpoints().size() == 3);
    assert(sub.get_triangles().size() == 1);
    assert(sub.get_quads().size() == 0);
    assert(sub.get_edges().size() == 0);

    auto vai = sub.get_vertex_atom_indices();
    assert(vai.size() == 3);
    for (const auto& atoms : vai) {
        assert(atoms.size() == 2);
        assert(colors[atoms[0]] != colors[atoms[1]]);
    }

    std::cout << "  3 vertices, 1 triangle - PASS\n";
}

void test_simplified_free_triangle() {
    std::cout << "Test: Simplified Free Triangle (2-1)\n";

    Points points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // color 1
        {0.5, 1.0, 0.0}                       // color 2
    };
    ColorLabels colors = {1, 1, 2};

    SimplifiedSubdivision sub(points, colors);
    sub.process_simplex({0, 1, 2});

    // 2-1 triangle: 2 cross-color edges → 2 vertices, 1 edge
    assert(sub.get_midpoints().size() == 2);
    assert(sub.get_triangles().size() == 0);
    assert(sub.get_quads().size() == 0);
    assert(sub.get_edges().size() == 1);

    std::cout << "  2 vertices, 1 edge - PASS\n";
}

void test_simplified_free_edge() {
    std::cout << "Test: Simplified Free Edge\n";

    Points points = {
        {0.0, 0.0, 0.0},  // color 1
        {1.0, 0.0, 0.0}   // color 2
    };
    ColorLabels colors = {1, 2};

    SimplifiedSubdivision sub(points, colors);
    sub.process_simplex({0, 1});

    // 1-1 edge: 1 cross-color pair → 1 vertex, no faces
    assert(sub.get_midpoints().size() == 1);
    assert(sub.get_triangles().size() == 0);
    assert(sub.get_quads().size() == 0);
    assert(sub.get_edges().size() == 0);

    std::cout << "  1 vertex - PASS\n";
}

void test_simplified_vertex_dedup() {
    std::cout << "Test: Simplified Vertex Deduplication\n";

    // Two tetrahedra sharing a multicolored face should share vertices
    Points points = {
        {0.0, 0.0, 0.0},   // 0: color 1
        {1.0, 0.0, 0.0},   // 1: color 1
        {0.5, 1.0, 0.0},   // 2: color 2
        {0.5, 0.5, 1.0},   // 3: color 2
        {0.5, 0.5, -1.0}   // 4: color 2
    };
    ColorLabels colors = {1, 1, 2, 2, 2};

    SimplifiedSubdivision sub(points, colors);
    // Tet 1: {0,1,2,3} — 2-2 → 4 vertices, 1 quad
    sub.process_tetrahedron({0, 1, 2, 3});
    assert(sub.get_midpoints().size() == 4);

    // Tet 2: {0,1,2,4} — 2-2, shares face {0,1,2} with tet 1
    // Cross-color pairs: (0,2), (0,4), (1,2), (1,4)
    // (0,2) and (1,2) already exist from tet 1 → only 2 new vertices
    sub.process_tetrahedron({0, 1, 2, 4});
    assert(sub.get_midpoints().size() == 6);  // 4 + 2 new

    // Total quads: 1 from each tet = 2
    assert(sub.get_quads().size() == 2);

    std::cout << "  6 vertices (4+2 with dedup), 2 quads - PASS\n";
}

void test_simplified_vs_barycentric_vertices() {
    std::cout << "Test: Simplified Vertices Match Barycentric Level-2\n";

    Points points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0}, {0.5, 0.5, 1.0}
    };
    ColorLabels colors = {1, 1, 2, 2};

    // Barycentric subdivision
    BarycentricSubdivision bary(points, colors);
    bary.process_tetrahedron({0, 1, 2, 3});
    auto bary_vai = bary.get_vertex_atom_indices();

    // Simplified subdivision
    SimplifiedSubdivision simp(points, colors);
    simp.process_tetrahedron({0, 1, 2, 3});
    auto simp_vai = simp.get_vertex_atom_indices();

    // Every simplified vertex (2 atoms) should exist in barycentric output
    for (const auto& simp_atoms : simp_vai) {
        assert(simp_atoms.size() == 2);
        bool found = false;
        for (const auto& bary_atoms : bary_vai) {
            if (bary_atoms == simp_atoms) {
                found = true;
                break;
            }
        }
        assert(found);
    }

    // Count 2-atom vertices in barycentric = number of simplified vertices
    size_t bary_level2 = 0;
    for (const auto& atoms : bary_vai) {
        if (atoms.size() == 2) bary_level2++;
    }
    assert(bary_level2 == simp_vai.size());

    // Verify that the actual 3D positions match
    auto bary_pts = bary.get_barycenters();
    auto simp_pts = simp.get_midpoints();

    for (size_t si = 0; si < simp_vai.size(); ++si) {
        for (size_t bi = 0; bi < bary_vai.size(); ++bi) {
            if (bary_vai[bi] == simp_vai[si]) {
                double dist = (bary_pts[bi] - simp_pts[si]).norm();
                assert(dist < 1e-10);
                break;
            }
        }
    }

    std::cout << "  All level-2 vertices match - PASS\n";
}

void test_simplified_rejects_3_colors() {
    std::cout << "Test: Simplified Rejects 3+ Colors\n";

    Points points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0}, {0.5, 0.5, 1.0}
    };
    ColorLabels colors = {1, 2, 3, 3};

    SimplifiedSubdivision sub(points, colors);

    bool caught = false;
    try {
        sub.process_tetrahedron({0, 1, 2, 3});
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    assert(caught);

    std::cout << "  Correctly rejects 3-color input - PASS\n";
}

void test_compute_simplified_surface_delaunay() {
    std::cout << "Test: compute_simplified_surface (unweighted Delaunay)\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0},
        {0.5, 0.5, 1.0},
        {2.0, 0.0, 0.0},
        {2.5, 1.0, 0.0}
    };
    ColorLabels colors = {1, 1, 1, 2, 2, 2};

    auto surface = compute_simplified_surface(points, colors, Radii{}, false);

    // Deterministic fixture: 3 generating tets yield 6 cross-color midpoints,
    // 2 triangles (3-1 tets) and 1 quad (2-2 tet).
    assert(surface.vertices.size() == 6);
    assert(surface.triangles.size() == 2);
    assert(surface.quads.size() == 1);
    assert(surface.edges.size() == 0);
    assert(surface.generating_simplices.generating_tetrahedra.size() == 3);
    assert(surface.weighted == false);
    assert(surface.alpha == false);

    // Every vertex should have exactly 2 atom indices
    for (const auto& atoms : surface.vertex_atom_indices) {
        assert(atoms.size() == 2);
        assert(colors[atoms[0]] != colors[atoms[1]]);
    }

    // The simplified surface must be strictly coarser than the barycentric one.
    auto bary_surface = compute_barycentric_subdivision_and_filtration(
        points, colors, Radii{}, false);
    auto& [bary_verts, bary_filt, bary_simp, bary_vai] = bary_surface;
    assert(bary_verts.size() == 17);
    assert(surface.vertices.size() < bary_verts.size());

    std::cout << "  PASS\n";
}

void test_compute_simplified_surface_weighted_alpha() {
    std::cout << "Test: compute_simplified_surface (weighted alpha)\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {2.0, 0.0, 0.0},
        {2.0, 1.0, 0.0}
    };
    ColorLabels colors = {1, 1, 1, 2, 2, 2};
    // Radii large enough that the alpha complex actually contains
    // multicolored simplices (0.3 left it empty and the test vacuous).
    Radii radii(points.size(), 1.0);

    auto surface = compute_simplified_surface(points, colors, radii, true);

    // Deterministic fixture: one 3-1 tet triangle plus two free bicolored
    // edges survive the alpha filtering at r=1.0.
    assert(surface.vertices.size() == 6);
    assert(surface.triangles.size() == 1);
    assert(surface.quads.size() == 0);
    assert(surface.edges.size() == 2);
    assert(surface.weighted == true);
    assert(surface.alpha == true);

    for (const auto& atoms : surface.vertex_atom_indices) {
        assert(atoms.size() == 2);
    }

    std::cout << "  PASS\n";
}

void test_quad_cyclic_order() {
    std::cout << "Test: Quad Vertices in Cyclic Order\n";

    // 2-2 tet: partition {a0,a1} vs {b0,b1}
    // Quad should be in cyclic order: a0b0 - a0b1 - a1b1 - a1b0
    // (consecutive vertices share one original atom)
    Points points = {
        {0.0, 0.0, 0.0},   // 0: color 1 (a0)
        {0.0, 3.0, 0.0},   // 1: color 1 (a1)
        {1.0, 0.0, 0.0},   // 2: color 2 (b0)
        {1.0, 1.0, 0.0}    // 3: color 2 (b1)
    };
    ColorLabels colors = {1, 1, 2, 2};

    SimplifiedSubdivision sub(points, colors);
    sub.process_tetrahedron({0, 1, 2, 3});

    assert(sub.get_quads().size() == 1);
    auto vai = sub.get_vertex_atom_indices();
    const auto& quad = sub.get_quads()[0];

    // Check that consecutive vertices in the quad share an original atom
    for (int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        const auto& atoms_i = vai[quad[i]];
        const auto& atoms_j = vai[quad[j]];
        bool shared = (atoms_i[0] == atoms_j[0] || atoms_i[0] == atoms_j[1] ||
                       atoms_i[1] == atoms_j[0] || atoms_i[1] == atoms_j[1]);
        assert(shared);
    }

    std::cout << "  Consecutive quad vertices share an atom - PASS\n";
}

// --- Helpers for random tests ---

static double triangle_area(const Point3D& p0, const Point3D& p1, const Point3D& p2) {
    return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

static std::pair<Points, ColorLabels> make_random_2color_cloud(int n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-5.0, 5.0);
    Points points(n);
    ColorLabels colors(n);
    for (int i = 0; i < n; ++i) {
        points[i] = Point3D(coord(rng), coord(rng), coord(rng));
        colors[i] = (i < n / 2) ? 0 : 1;
    }
    return {points, colors};
}

void test_random_vertices_equal_multicolored_edges() {
    std::cout << "Test: Random Cloud — Vertices = Multicolored Delaunay Edges\n";

    auto [points, colors] = make_random_2color_cloud(60, /*seed=*/42);

    // Get all multicolored tetrahedra from the Delaunay triangulation
    InterfaceGenerator gen;
    auto tets = gen.collect_multicolored_simplices(points, colors, Radii{}, false).generating_tetrahedra;

    // Collect unique multicolored edges
    std::set<std::array<int, 2>> mc_edges;
    for (const auto& tet : tets) {
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                if (colors[tet[i]] != colors[tet[j]]) {
                    int a = std::min(tet[i], tet[j]);
                    int b = std::max(tet[i], tet[j]);
                    mc_edges.insert({a, b});
                }
            }
        }
    }

    // Compute simplified surface
    auto surface = compute_simplified_surface(points, colors, Radii{}, false);

    std::cout << "  Multicolored Delaunay edges: " << mc_edges.size() << "\n";
    std::cout << "  Simplified vertices: " << surface.vertices.size() << "\n";
    assert(surface.vertices.size() == mc_edges.size());

    std::cout << "  PASS\n";
}

void test_random_area_equality() {
    std::cout << "Test: Random Cloud — Simplified Area = Barycentric Area\n";

    auto [points, colors] = make_random_2color_cloud(60, /*seed=*/42);

    // Barycentric surface
    auto [bary_verts, bary_filt, bary_simp, bary_vai] =
        compute_barycentric_subdivision_and_filtration(points, colors, Radii{}, false);

    double bary_area = 0.0;
    for (const auto& [simplex, val] : bary_filt) {
        if (simplex.size() == 3) {
            bary_area += triangle_area(
                bary_verts[simplex[0]], bary_verts[simplex[1]], bary_verts[simplex[2]]);
        }
    }

    // Simplified surface
    auto surface = compute_simplified_surface(points, colors, Radii{}, false);

    double simp_area = 0.0;
    for (const auto& tri : surface.triangles) {
        simp_area += triangle_area(
            surface.vertices[tri[0]], surface.vertices[tri[1]], surface.vertices[tri[2]]);
    }
    for (const auto& quad : surface.quads) {
        // Quad in cyclic order → split along diagonal 0-2
        simp_area += triangle_area(
            surface.vertices[quad[0]], surface.vertices[quad[1]], surface.vertices[quad[2]]);
        simp_area += triangle_area(
            surface.vertices[quad[0]], surface.vertices[quad[2]], surface.vertices[quad[3]]);
    }

    std::cout << "  Barycentric area:  " << bary_area << "\n";
    std::cout << "  Simplified area:   " << simp_area << "\n";
    std::cout << "  Relative diff:     " << std::abs(bary_area - simp_area) / bary_area << "\n";
    assert(std::abs(bary_area - simp_area) < 1e-10 * bary_area);

    std::cout << "  PASS\n";
}

void test_random_mesh_manifold() {
    std::cout << "Test: Random Cloud — Mesh is Manifold (no duplicate/non-manifold edges)\n";

    auto [points, colors] = make_random_2color_cloud(60, /*seed=*/42);
    auto surface = compute_simplified_surface(points, colors, Radii{}, false);

    using Edge = std::array<int32_t, 2>;
    auto make_edge = [](int32_t a, int32_t b) -> Edge {
        return {std::min(a, b), std::max(a, b)};
    };

    // Count how many faces each edge borders
    std::map<Edge, int> edge_count;

    for (const auto& tri : surface.triangles) {
        edge_count[make_edge(tri[0], tri[1])]++;
        edge_count[make_edge(tri[1], tri[2])]++;
        edge_count[make_edge(tri[0], tri[2])]++;
    }
    for (const auto& quad : surface.quads) {
        edge_count[make_edge(quad[0], quad[1])]++;
        edge_count[make_edge(quad[1], quad[2])]++;
        edge_count[make_edge(quad[2], quad[3])]++;
        edge_count[make_edge(quad[3], quad[0])]++;
    }

    int boundary = 0, interior = 0, non_manifold = 0;
    for (const auto& [edge, count] : edge_count) {
        if (count == 1) boundary++;
        else if (count == 2) interior++;
        else non_manifold++;
    }

    std::cout << "  Total mesh edges: " << edge_count.size() << "\n";
    std::cout << "  Interior (count=2): " << interior << "\n";
    std::cout << "  Boundary (count=1): " << boundary << "\n";
    std::cout << "  Non-manifold (count>2): " << non_manifold << "\n";
    assert(non_manifold == 0);

    // Check no duplicate faces
    std::set<std::vector<int32_t>> face_set;
    for (const auto& tri : surface.triangles) {
        std::vector<int32_t> key(tri.begin(), tri.end());
        std::sort(key.begin(), key.end());
        assert(face_set.insert(key).second);  // insert must succeed (no duplicate)
    }
    for (const auto& quad : surface.quads) {
        std::vector<int32_t> key(quad.begin(), quad.end());
        std::sort(key.begin(), key.end());
        assert(face_set.insert(key).second);
    }

    std::cout << "  Faces: " << face_set.size()
              << " (" << surface.triangles.size() << " tri + "
              << surface.quads.size() << " quad), all unique\n";

    // The fixture is deterministic (fixed seed), so counts are exact.
    int V = static_cast<int>(surface.vertices.size());
    int E = static_cast<int>(edge_count.size());
    int F = static_cast<int>(face_set.size());
    assert(V == 200);
    assert(E == 475);
    assert(F == 266);
    assert(V - E + F == -9);

    std::cout << "  PASS\n";
}

int main() {
    std::cout << "Running Simplified Subdivision Tests\n";
    std::cout << "====================================\n\n";

    try {
        test_simplified_22_tetrahedron();
        test_simplified_31_tetrahedron();
        test_simplified_free_triangle();
        test_simplified_free_edge();
        test_simplified_vertex_dedup();
        test_simplified_vs_barycentric_vertices();
        test_simplified_rejects_3_colors();
        test_compute_simplified_surface_delaunay();
        test_compute_simplified_surface_weighted_alpha();
        test_quad_cyclic_order();
        test_random_vertices_equal_multicolored_edges();
        test_random_area_equality();
        test_random_mesh_manifold();

        std::cout << "\nAll simplified subdivision tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
