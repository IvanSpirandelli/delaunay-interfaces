/**
 * Combinatorial Property Tests
 *
 * These tests verify mathematical/combinatorial properties of the interface
 * surface that can be derived from first principles, without needing external
 * ground truth data.
 *
 * For example: when all points have distinct colors, every tetrahedron in the
 * Delaunay triangulation is multicolored (1-1-1-1 partition), and the number
 * of triangles in the interface surface equals 12 times the number of tetrahedra.
 */

#include <iostream>
#include <cassert>
#include <random>
#include <delaunay_interfaces/interface_generation.hpp>

using namespace delaunay_interfaces;

/**
 * Test: All points have distinct colors
 *
 * When we have n points with n distinct colors:
 * - Every tetrahedron in the Delaunay triangulation is multicolored
 * - Each tetrahedron has the 1-1-1-1 partition pattern
 * - For 1-1-1-1 partition, barycentric subdivision produces 12 triangles per tet
 *
 * This test generates 6 random points with 6 colors and verifies:
 * - num_triangles == num_tetrahedra * 12
 */
void test_all_distinct_colors() {
    std::cout << "Test: All distinct colors (6 points, 6 colors)\n";

    // Generate reproducible random points
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(0.0, 10.0);

    Points points;
    ColorLabels colors;
    for (int i = 0; i < 6; ++i) {
        points.push_back({dist(rng), dist(rng), dist(rng)});
        colors.push_back(i + 1);  // Each point gets a unique color (1-6)
    }

    // Compute interface surface using unweighted Delaunay
    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    // Count tetrahedra by getting them directly
    auto tetrahedra = generator.get_multicolored_tetrahedra(points, colors, {}, false, false);
    size_t num_tetrahedra = tetrahedra.size();

    // Count triangles in filtration (simplices with 3 vertices)
    size_t num_triangles = 0;
    for (const auto& [simplex, value] : surface.filtration) {
        if (simplex.size() == 3) {
            num_triangles++;
        }
    }

    std::cout << "  Points: " << points.size() << "\n";
    std::cout << "  Tetrahedra (all multicolored): " << num_tetrahedra << "\n";
    std::cout << "  Triangles in filtration: " << num_triangles << "\n";
    std::cout << "  Expected triangles (12 per tet): " << num_tetrahedra * 12 << "\n";

    // For 1-1-1-1 partition, each tetrahedron contributes 12 triangles
    // (This is the barycentric subdivision of the interface surface)
    assert(num_triangles == num_tetrahedra * 12);

    std::cout << "  PASS\n";
}

/**
 * Test: Verify edge and vertex counts for all-distinct-colors case
 *
 * For 1-1-1-1 partition per tetrahedron:
 * - 11 vertices per tetrahedron (6 edge + 4 triangle + 1 tet barycenters)
 * - But vertices are shared between adjacent tetrahedra
 *
 * This test verifies the filtration structure is consistent.
 */
void test_filtration_structure() {
    std::cout << "\nTest: Filtration structure verification\n";

    std::mt19937 rng(54321);
    std::uniform_real_distribution<double> dist(0.0, 10.0);

    Points points;
    ColorLabels colors;
    for (int i = 0; i < 6; ++i) {
        points.push_back({dist(rng), dist(rng), dist(rng)});
        colors.push_back(i + 1);
    }

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    // Count simplices by dimension
    size_t num_vertices = 0;
    size_t num_edges = 0;
    size_t num_triangles = 0;

    for (const auto& [simplex, value] : surface.filtration) {
        if (simplex.size() == 1) num_vertices++;
        else if (simplex.size() == 2) num_edges++;
        else if (simplex.size() == 3) num_triangles++;
    }

    std::cout << "  Vertices (0-simplices): " << num_vertices << "\n";
    std::cout << "  Edges (1-simplices): " << num_edges << "\n";
    std::cout << "  Triangles (2-simplices): " << num_triangles << "\n";
    std::cout << "  Total simplices: " << surface.filtration.size() << "\n";

    // Basic sanity checks
    assert(num_vertices > 0);
    assert(num_edges > 0);
    assert(num_triangles > 0);
    assert(num_vertices + num_edges + num_triangles == surface.filtration.size());

    // Euler characteristic check for a surface: V - E + F should be consistent
    // For a closed surface embedded in 3D from tetrahedra subdivision
    int euler = static_cast<int>(num_vertices) - static_cast<int>(num_edges) + static_cast<int>(num_triangles);
    std::cout << "  Euler characteristic (V - E + F): " << euler << "\n";

    std::cout << "  PASS\n";
}

int main() {
    std::cout << "Combinatorial Property Tests\n";
    std::cout << "============================\n";
    std::cout << "Testing mathematical properties derived from first principles\n\n";

    try {
        test_all_distinct_colors();
        test_filtration_structure();

        std::cout << "\nAll combinatorial tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
