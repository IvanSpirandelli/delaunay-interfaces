/**
 * Custom barycenter verification test
 *
 * Verifies barycenter computation for a simple unit tetrahedron with 4 distinct colors.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <delaunay_interfaces/interface_generation.hpp>

using namespace delaunay_interfaces;

bool approx_equal(double a, double b, double tol = 1e-6) {
    return std::abs(a - b) < tol;
}

bool point_approx_equal(const Point3D& p, double x, double y, double z, double tol = 1e-6) {
    return approx_equal(p.x(), x, tol) && approx_equal(p.y(), y, tol) && approx_equal(p.z(), z, tol);
}

void test_3_1_partition() {
    std::cout << "\nTest: 3-1 partition (equilateral triangle + apex)\n";
    std::cout << "=================================================\n";

    // Equilateral triangle in x-z plane at y = -1/3, apex at y = 1
    // User's original specification
    double y_triangle = -1.0 / 3.0;
    double y_apex = 1.0;

    // Equilateral triangle vertices (radius 1 from origin in x-z plane)
    double r = 1.0;
    Points points = {
        {r, y_triangle, 0.0},                          // p0
        {-r/2.0, y_triangle, r * std::sqrt(3.0)/2.0},  // p1
        {-r/2.0, y_triangle, -r * std::sqrt(3.0)/2.0}, // p2
        {0.0, y_apex, 0.0}                              // p3 (apex)
    };
    ColorLabels colors = {1, 1, 1, 2};  // 3-1 partition

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    std::cout << "Input points:\n";
    for (size_t i = 0; i < points.size(); ++i) {
        std::cout << "  p" << i << " = [" << std::fixed << std::setprecision(4)
                  << points[i].x() << ", " << points[i].y() << ", " << points[i].z()
                  << "] color=" << colors[i] << "\n";
    }

    std::cout << "\nComputed vertices (" << surface.vertices.size() << " total):\n";
    std::cout << std::fixed << std::setprecision(6);

    double expected_edge_y = (y_triangle + y_apex) / 2.0;
    std::cout << "\nExpected edge barycenter y = (" << y_triangle << " + " << y_apex << ")/2 = " << expected_edge_y << "\n";

    for (const auto& v : surface.vertices) {
        std::cout << "  [" << v.x() << ", " << v.y() << ", " << v.z() << "]\n";
    }

    // Count triangles
    size_t num_triangles = 0;
    for (const auto& [simplex, value] : surface.filtration) {
        if (simplex.size() == 3) num_triangles++;
    }
    std::cout << "\nNumber of triangles: " << num_triangles << "\n";

    // For 3-1 partition, we expect 6 triangles
    std::cout << "Expected triangles for 3-1 partition: 6\n";
    if (num_triangles == 6) {
        std::cout << "PASS (triangle count)\n";
    } else {
        std::cout << "FAIL: Expected 6 triangles!\n";
    }
}

void test_3_1_partition_y_zero() {
    std::cout << "\nTest: 3-1 partition with interface at y=0\n";
    std::cout << "==========================================\n";

    // For interface at y=0: triangle at y=-1, apex at y=1
    // Edge barycenters: (-1 + 1)/2 = 0
    double y_triangle = -1.0;
    double y_apex = 1.0;

    double r = 1.0;
    Points points = {
        {r, y_triangle, 0.0},
        {-r/2.0, y_triangle, r * std::sqrt(3.0)/2.0},
        {-r/2.0, y_triangle, -r * std::sqrt(3.0)/2.0},
        {0.0, y_apex, 0.0}
    };
    ColorLabels colors = {1, 1, 1, 2};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    std::cout << "Input points (triangle at y=-1, apex at y=1):\n";
    for (size_t i = 0; i < points.size(); ++i) {
        std::cout << "  p" << i << " = [" << std::fixed << std::setprecision(4)
                  << points[i].x() << ", " << points[i].y() << ", " << points[i].z()
                  << "] color=" << colors[i] << "\n";
    }

    std::cout << "\nComputed vertices (" << surface.vertices.size() << " total):\n";

    bool all_at_y_zero = true;
    for (const auto& v : surface.vertices) {
        std::cout << "  [" << std::fixed << std::setprecision(6)
                  << v.x() << ", " << v.y() << ", " << v.z() << "]";
        if (!approx_equal(v.y(), 0.0)) {
            std::cout << " <-- y != 0!";
            all_at_y_zero = false;
        }
        std::cout << "\n";
    }

    size_t num_triangles = 0;
    for (const auto& [simplex, value] : surface.filtration) {
        if (simplex.size() == 3) num_triangles++;
    }

    std::cout << "\nVerification:\n";
    std::cout << "  All vertices at y=0: " << (all_at_y_zero ? "YES" : "NO") << "\n";
    std::cout << "  Number of triangles: " << num_triangles << " (expected 6)\n";

    if (all_at_y_zero && num_triangles == 6) {
        std::cout << "  PASS\n";
    } else {
        if (!all_at_y_zero) std::cout << "  FAIL: Interface should lie at y=0!\n";
        if (num_triangles != 6) std::cout << "  FAIL: Expected 6 triangles!\n";
        exit(1);
    }
}

void test_2_2_partition_y_plane() {
    std::cout << "\nTest: 2-2 partition with interface at y=0\n";
    std::cout << "==========================================\n";

    // Two points at y=-1 (color 1), two at y=1 (color 2)
    // Interface should lie entirely at y=0
    Points points = {
        {-1.0, -1.0, 0.0},  // color 1
        {1.0, -1.0, 0.0},   // color 1
        {0.0, 1.0, -1.0},   // color 2
        {0.0, 1.0, 1.0}     // color 2
    };
    ColorLabels colors = {1, 1, 2, 2};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    std::cout << "Input points:\n";
    for (size_t i = 0; i < points.size(); ++i) {
        std::cout << "  p" << i << " = [" << points[i].x() << ", " << points[i].y() << ", " << points[i].z()
                  << "] color=" << colors[i] << "\n";
    }

    std::cout << "\nComputed vertices (" << surface.vertices.size() << " total):\n";
    std::cout << std::fixed << std::setprecision(6);

    bool all_at_y_zero = true;
    for (const auto& v : surface.vertices) {
        std::cout << "  [" << v.x() << ", " << v.y() << ", " << v.z() << "]";
        if (!approx_equal(v.y(), 0.0)) {
            std::cout << " <-- y != 0!";
            all_at_y_zero = false;
        }
        std::cout << "\n";
    }

    // Count triangles in filtration
    size_t num_triangles = 0;
    for (const auto& [simplex, value] : surface.filtration) {
        if (simplex.size() == 3) {
            num_triangles++;
        }
    }

    std::cout << "\nVerification:\n";
    std::cout << "  All vertices at y=0: " << (all_at_y_zero ? "YES" : "NO") << "\n";
    std::cout << "  Number of triangles: " << num_triangles << "\n";

    if (all_at_y_zero && num_triangles == 8) {
        std::cout << "  PASS\n";
    } else {
        if (!all_at_y_zero) std::cout << "  FAIL: Interface should lie at y=0!\n";
        if (num_triangles != 8) std::cout << "  FAIL: Expected 8 triangles, got " << num_triangles << "!\n";
        exit(1);
    }
}

int main() {
    std::cout << "Custom Barycenter Verification Test\n";
    std::cout << "====================================\n\n";

    // Test 3-1 partitions
    test_3_1_partition();
    test_3_1_partition_y_zero();

    // Test 2-2 partition
    test_2_2_partition_y_plane();

    std::cout << "\n";

    // Unit tetrahedron with 4 distinct colors (1-1-1-1 partition)
    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    ColorLabels colors = {1, 2, 3, 4};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

    std::cout << "Input points:\n";
    for (size_t i = 0; i < points.size(); ++i) {
        std::cout << "  p" << i << " = [" << points[i][0] << ", " << points[i][1] << ", " << points[i][2]
                  << "] color=" << colors[i] << "\n";
    }

    std::cout << "\nComputed vertices (" << surface.vertices.size() << " total):\n";
    std::cout << std::fixed << std::setprecision(6);

    // Sort vertices for consistent output
    std::vector<std::array<double, 3>> sorted_verts;
    for (const auto& v : surface.vertices) {
        sorted_verts.push_back({v.x(), v.y(), v.z()});
    }
    std::sort(sorted_verts.begin(), sorted_verts.end());

    for (const auto& v : sorted_verts) {
        std::cout << "  [" << v[0] << ", " << v[1] << ", " << v[2] << "]\n";
    }

    // Expected barycenters for 1-1-1-1 partition:
    // 6 edge barycenters
    std::cout << "\nExpected edge barycenters (6):\n";
    std::cout << "  [0.5, 0.0, 0.0] - edge p0-p1\n";
    std::cout << "  [0.0, 0.5, 0.0] - edge p0-p2\n";
    std::cout << "  [0.0, 0.0, 0.5] - edge p0-p3\n";
    std::cout << "  [0.5, 0.5, 0.0] - edge p1-p2\n";
    std::cout << "  [0.5, 0.0, 0.5] - edge p1-p3\n";
    std::cout << "  [0.0, 0.5, 0.5] - edge p2-p3\n";

    // 4 triangle barycenters
    double third = 1.0 / 3.0;
    std::cout << "\nExpected triangle barycenters (4):\n";
    std::cout << "  [" << third << ", " << third << ", 0.0] - face p0-p1-p2\n";
    std::cout << "  [" << third << ", 0.0, " << third << "] - face p0-p1-p3\n";
    std::cout << "  [0.0, " << third << ", " << third << "] - face p0-p2-p3\n";
    std::cout << "  [" << third << ", " << third << ", " << third << "] - face p1-p2-p3\n";

    // 1 tetrahedron barycenter
    std::cout << "\nExpected tetrahedron barycenter (1):\n";
    std::cout << "  [0.25, 0.25, 0.25]\n";

    // Verify counts
    std::cout << "\nVerification:\n";
    std::cout << "  Expected vertices: 11 (6 edge + 4 triangle + 1 tet)\n";
    std::cout << "  Computed vertices: " << surface.vertices.size() << "\n";

    if (surface.vertices.size() != 11) {
        std::cout << "  FAIL: Vertex count mismatch!\n";
        return 1;
    }

    // Check specific expected points
    bool all_found = true;
    std::vector<std::array<double, 3>> expected = {
        {0.5, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 0.5},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5},
        {third, third, 0.0},
        {third, 0.0, third},
        {0.0, third, third},
        {third, third, third},
        {0.25, 0.25, 0.25}
    };

    for (const auto& exp : expected) {
        bool found = false;
        for (const auto& v : surface.vertices) {
            if (point_approx_equal(v, exp[0], exp[1], exp[2])) {
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "  MISSING: [" << exp[0] << ", " << exp[1] << ", " << exp[2] << "]\n";
            all_found = false;
        }
    }

    if (all_found) {
        std::cout << "  All expected barycenters found!\n";
        std::cout << "\nPASS\n";
        return 0;
    } else {
        std::cout << "\nFAIL: Some barycenters missing!\n";
        return 1;
    }
}
