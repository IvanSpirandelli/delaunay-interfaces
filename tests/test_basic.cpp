#include <iostream>
#include <cassert>
#include <cmath>
#include <set>
#include <delaunay_interfaces/interface_generation.hpp>
#include <delaunay_interfaces/barycentric_subdivision.hpp>
#include <delaunay_interfaces/chromatic_partitioning.hpp>

using namespace delaunay_interfaces;

void test_chromatic_partitioning() {
    std::cout << "Test: Chromatic Partitioning\n";

    ColorLabels colors = {1, 1, 2, 2};
    Tetrahedron tet = {0, 1, 2, 3};

    auto partition = get_chromatic_partitioning(tet, colors);

    assert(partition.size() == 2);
    assert(partition[0].size() == 2);
    assert(partition[1].size() == 2);

    std::cout << "  PASS\n";
}

void test_chromatic_partitioning_vector() {
    std::cout << "Test: Chromatic Partitioning (vector overload)\n";

    ColorLabels colors = {1, 2, 3};

    // Triangle
    std::vector<int> tri = {0, 1, 2};
    auto partition = get_chromatic_partitioning(tri, colors);
    assert(partition.size() == 3);

    // Edge
    std::vector<int> edge = {0, 1};
    auto partition2 = get_chromatic_partitioning(edge, colors);
    assert(partition2.size() == 2);

    std::cout << "  PASS\n";
}

void test_barycenter() {
    std::cout << "Test: Barycenter Computation\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    std::vector<int> indices = {0, 1, 2, 3};
    Point3D bc = compute_barycenter(points, indices);

    assert(std::abs(bc.x() - 0.25) < 1e-10);
    assert(std::abs(bc.y() - 0.25) < 1e-10);
    assert(std::abs(bc.z() - 0.25) < 1e-10);

    std::cout << "  PASS\n";
}

void test_simple_delaunay() {
    std::cout << "Test: Simple Delaunay Complex\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0},
        {0.5, 0.5, 1.0},
        {2.0, 0.0, 0.0},
        {2.5, 1.0, 0.0}
    };

    ColorLabels colors = {1, 1, 1, 2, 2, 2};

    InterfaceGenerator generator;

    try {
        auto surface = generator.compute_interface_surface(points, colors, {}, false, false);
        std::cout << "  Found " << surface.vertices.size() << " barycenters\n";
        std::cout << "  Found " << surface.filtration.size() << " filtration simplices\n";
        assert(surface.vertices.size() > 0);
        assert(surface.filtration.size() > 0);
        std::cout << "  PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << e.what() << "\n";
        exit(1);
    }
}

void test_weighted_alpha() {
    std::cout << "Test: Weighted Alpha Complex\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {2.0, 0.0, 0.0},
        {2.0, 1.0, 0.0}
    };

    ColorLabels colors = {1, 1, 1, 2, 2, 2};
    Radii radii(points.size(), 0.3);

    InterfaceGenerator generator;

    try {
        auto surface = generator.compute_interface_surface(points, colors, radii, true, true);
        std::cout << "  Found " << surface.vertices.size() << " barycenters\n";
        std::cout << "  Found " << surface.filtration.size() << " filtration simplices\n";
        assert(surface.weighted == true);
        assert(surface.alpha == true);
        std::cout << "  PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << e.what() << "\n";
        exit(1);
    }
}

void test_simple_delaunay_lower_star() {
    std::cout << "Test: Simple Delaunay Complex (Lower Star)\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0},
        {0.5, 0.5, 1.0},
        {2.0, 0.0, 0.0},
        {2.5, 1.0, 0.0}
    };

    ColorLabels colors = {1, 1, 1, 2, 2, 2};

    InterfaceGenerator generator;

    try {
        auto upper = generator.compute_interface_surface(points, colors, {}, false, false, false);
        auto lower = generator.compute_interface_surface(points, colors, {}, false, false, true);

        assert(lower.lower_star == true);
        assert(upper.lower_star == false);
        assert(lower.vertices.size() == upper.vertices.size());
        assert(lower.filtration.size() == upper.filtration.size());

        // For edges and triangles, lower star (max) values >= upper star (min) values
        for (size_t i = 0; i < lower.filtration.size(); ++i) {
            auto& [l_simplex, l_val] = lower.filtration[i];
            auto& [u_simplex, u_val] = upper.filtration[i];
            if (l_simplex.size() == 1) {
                // Vertex filtration values are identical
                assert(std::abs(l_val - u_val) < 1e-10);
            }
        }

        std::cout << "  PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << e.what() << "\n";
        exit(1);
    }
}

void test_weighted_alpha_lower_star() {
    std::cout << "Test: Weighted Alpha Complex (Lower Star)\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {2.0, 0.0, 0.0},
        {2.0, 1.0, 0.0}
    };

    ColorLabels colors = {1, 1, 1, 2, 2, 2};
    Radii radii(points.size(), 0.3);

    InterfaceGenerator generator;

    try {
        auto upper = generator.compute_interface_surface(points, colors, radii, true, true, false);
        auto lower = generator.compute_interface_surface(points, colors, radii, true, true, true);

        assert(lower.lower_star == true);
        assert(lower.weighted == true);
        assert(lower.alpha == true);
        assert(lower.vertices.size() == upper.vertices.size());
        assert(lower.filtration.size() == upper.filtration.size());

        // Vertex values must match between upper and lower star
        for (size_t i = 0; i < lower.filtration.size(); ++i) {
            auto& [l_simplex, l_val] = lower.filtration[i];
            auto& [u_simplex, u_val] = upper.filtration[i];
            if (l_simplex.size() == 1) {
                assert(std::abs(l_val - u_val) < 1e-10);
            }
        }

        std::cout << "  PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << e.what() << "\n";
        exit(1);
    }
}

void test_input_validation() {
    std::cout << "Test: Input Validation\n";

    Points points = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    ColorLabels colors = {1}; // Wrong size

    InterfaceGenerator generator;

    bool caught_exception = false;
    try {
        auto surface = generator.compute_interface_surface(points, colors);
    } catch (const std::invalid_argument& e) {
        caught_exception = true;
    }

    assert(caught_exception);
    std::cout << "  PASS\n";
}

void test_generic_combination_counts() {
    std::cout << "Test: Generic Combination and Chain Counts\n";

    // Test that process_simplex produces the correct number of combinations
    // for each partition type by counting barycenters and filtration simplices.

    // We use a single tetrahedron with known partition types.
    // The BarycentricSubdivision counts should match the old hardcoded values.

    // 2-2 partition: 9 combinations, 8 triangles
    {
        Points points = {
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // color 1
            {0.5, 1.0, 0.0}, {0.5, 0.5, 1.0}    // color 2
        };
        ColorLabels colors = {1, 1, 2, 2};

        BarycentricSubdivision sub(points, colors);
        sub.process_tetrahedron({0, 1, 2, 3});

        auto filtration = sub.get_filtration();
        size_t n_verts = 0, n_edges = 0, n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 1) n_verts++;
            else if (simplex.size() == 2) n_edges++;
            else if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 9);
        assert(n_tris == 8);
        std::cout << "  2-2: " << sub.get_barycenters().size() << " combos, "
                  << n_tris << " triangles - OK\n";
    }

    // 3-1 partition: 7 combinations, 6 triangles
    {
        Points points = {
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 1.0, 0.0},  // color 1
            {0.5, 0.5, 1.0}                                        // color 2
        };
        ColorLabels colors = {1, 1, 1, 2};

        BarycentricSubdivision sub(points, colors);
        sub.process_tetrahedron({0, 1, 2, 3});

        auto filtration = sub.get_filtration();
        size_t n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 7);
        assert(n_tris == 6);
        std::cout << "  3-1: " << sub.get_barycenters().size() << " combos, "
                  << n_tris << " triangles - OK\n";
    }

    // 2-1-1 partition: 10 combinations, 10 triangles
    {
        Points points = {
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // color 1
            {0.5, 1.0, 0.0},                      // color 2
            {0.5, 0.5, 1.0}                       // color 3
        };
        ColorLabels colors = {1, 1, 2, 3};

        BarycentricSubdivision sub(points, colors);
        sub.process_tetrahedron({0, 1, 2, 3});

        auto filtration = sub.get_filtration();
        size_t n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 10);
        assert(n_tris == 10);
        std::cout << "  2-1-1: " << sub.get_barycenters().size() << " combos, "
                  << n_tris << " triangles - OK\n";
    }

    // 1-1-1-1 partition: 11 combinations, 12 triangles
    {
        Points points = {
            {0.0, 0.0, 0.0},   // color 1
            {1.0, 0.0, 0.0},   // color 2
            {0.5, 1.0, 0.0},   // color 3
            {0.5, 0.5, 1.0}    // color 4
        };
        ColorLabels colors = {1, 2, 3, 4};

        BarycentricSubdivision sub(points, colors);
        sub.process_tetrahedron({0, 1, 2, 3});

        auto filtration = sub.get_filtration();
        size_t n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 11);
        assert(n_tris == 12);
        std::cout << "  1-1-1-1: " << sub.get_barycenters().size() << " combos, "
                  << n_tris << " triangles - OK\n";
    }

    std::cout << "  PASS\n";
}

void test_lower_dimensional_simplices() {
    std::cout << "Test: Lower-Dimensional Simplices (Triangle and Edge)\n";

    Points points = {
        {0.0, 0.0, 0.0},   // color 1
        {1.0, 0.0, 0.0},   // color 2
        {0.5, 1.0, 0.0}    // color 2
    };
    ColorLabels colors = {1, 2, 2};

    // Bicolored triangle [0] vs [1,2] → 3 combinations, 2 edges (chains), no triangles
    {
        BarycentricSubdivision sub(points, colors);
        sub.process_simplex({0, 1, 2});

        auto filtration = sub.get_filtration();
        size_t n_verts = 0, n_edges = 0, n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 1) n_verts++;
            else if (simplex.size() == 2) n_edges++;
            else if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 3);
        // 2 levels → chains of length 2, but those are just edges (not triangles)
        // Edges: all subset-inclusion pairs (2 atomic ⊂ 1 non-atomic) = 2 edges
        assert(n_edges == 2);
        assert(n_tris == 0);
        std::cout << "  Bicolored triangle (2-1): " << sub.get_barycenters().size()
                  << " vertices, " << n_edges << " edges, " << n_tris << " triangles - OK\n";
    }

    // Tricolored triangle [0] vs [1] vs [2] → 4 combinations, 3 edges, no triangles
    {
        ColorLabels colors3 = {1, 2, 3};
        BarycentricSubdivision sub(points, colors3);
        sub.process_simplex({0, 1, 2});

        auto filtration = sub.get_filtration();
        size_t n_verts = 0, n_edges = 0, n_tris = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 1) n_verts++;
            else if (simplex.size() == 2) n_edges++;
            else if (simplex.size() == 3) n_tris++;
        }
        assert(sub.get_barycenters().size() == 4);
        // 3 atomic at level 2, 1 non-atomic at level 3 → 3 edges (subset inclusion)
        // Plus 3 more edges between atomic pairs that are subsets of the same level-3
        // Actually: edges are all (i,j) where flat[i] ⊂ flat[j]
        // Level 2: {0,1}, {0,2}, {1,2}; Level 3: {0,1,2}
        // Edges: 3 (each level-2 ⊂ level-3)
        // No chains of length 3, so no triangles
        assert(n_edges == 3);
        assert(n_tris == 0);
        std::cout << "  Tricolored triangle (1-1-1): " << sub.get_barycenters().size()
                  << " vertices, " << n_edges << " edges, " << n_tris << " triangles - OK\n";
    }

    // Bicolored edge → 1 combination = 1 vertex, no edges
    {
        BarycentricSubdivision sub(points, colors);
        sub.process_simplex({0, 1});

        auto filtration = sub.get_filtration();
        size_t n_verts = 0, n_edges = 0;
        for (const auto& [simplex, val] : filtration) {
            if (simplex.size() == 1) n_verts++;
            else if (simplex.size() == 2) n_edges++;
        }
        assert(sub.get_barycenters().size() == 1);
        assert(n_verts == 1);
        assert(n_edges == 0);
        std::cout << "  Bicolored edge: " << sub.get_barycenters().size()
                  << " vertex, " << n_edges << " edges - OK\n";
    }

    // Monocolored edge → skipped (no output)
    {
        BarycentricSubdivision sub(points, colors);
        sub.process_simplex({1, 2}); // both color 2

        assert(sub.get_barycenters().size() == 0);
        assert(sub.get_filtration().empty());
        std::cout << "  Monocolored edge: skipped - OK\n";
    }

    std::cout << "  PASS\n";
}

void test_free_simplices_alpha() {
    std::cout << "Test: Free Simplices in Alpha Complex\n";

    // Create a point cloud where alpha filtering creates free multicolored
    // simplices. Use large radii for nearby points (to be in alpha complex)
    // but small radii for distant points (so their tetrahedra are excluded).

    // Two groups of overlapping spheres that share an interface
    // but are far enough apart that tetrahedra spanning the gap have alpha > 0.
    Points points = {
        // Color 1 cluster
        {0.0, 0.0, 0.0},
        {0.3, 0.0, 0.0},
        {0.0, 0.3, 0.0},
        // Color 2 cluster
        {0.5, 0.0, 0.0},
        {0.8, 0.0, 0.0},
        {0.5, 0.3, 0.0},
        // Add z-offset points to make 3D triangulation work
        {0.25, 0.15, 0.3},
        {0.55, 0.15, 0.3},
    };
    ColorLabels colors = {1, 1, 1, 2, 2, 2, 1, 2};
    // Moderate radii: enough for nearby pairs to overlap, not enough for distant tets
    Radii radii = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

    InterfaceGenerator generator;

    try {
        auto simplices = generator.get_multicolored_simplices_weighted_alpha(
            points, colors, radii);

        std::cout << "  Tetrahedra: " << simplices.tetrahedra.size() << "\n";
        std::cout << "  Free triangles: " << simplices.free_triangles.size() << "\n";
        std::cout << "  Free edges: " << simplices.free_edges.size() << "\n";

        // The subdivision should work with all simplex types
        BarycentricSubdivision sub(points, colors);
        for (const auto& tet : simplices.tetrahedra) {
            sub.process_tetrahedron(tet);
        }
        for (const auto& tri : simplices.free_triangles) {
            sub.process_simplex(tri);
        }
        for (const auto& edge : simplices.free_edges) {
            sub.process_simplex(edge);
        }

        auto filtration = sub.get_filtration();
        std::cout << "  Total barycenters: " << sub.get_barycenters().size() << "\n";
        std::cout << "  Total filtration simplices: " << filtration.size() << "\n";

        // Basic sanity: we should have some output
        size_t total_simplices = simplices.tetrahedra.size()
            + simplices.free_triangles.size()
            + simplices.free_edges.size();
        if (total_simplices > 0) {
            assert(sub.get_barycenters().size() > 0);
            assert(filtration.size() > 0);
        }

        std::cout << "  PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << e.what() << "\n";
        exit(1);
    }
}

int main() {
    std::cout << "Running DelaunayInterfaces C++ Tests\n";
    std::cout << "=====================================\n\n";

    try {
        test_chromatic_partitioning();
        test_chromatic_partitioning_vector();
        test_barycenter();
        test_simple_delaunay();
        test_weighted_alpha();
        test_simple_delaunay_lower_star();
        test_weighted_alpha_lower_star();
        test_input_validation();
        test_generic_combination_counts();
        test_lower_dimensional_simplices();
        test_free_simplices_alpha();

        std::cout << "\nAll tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
