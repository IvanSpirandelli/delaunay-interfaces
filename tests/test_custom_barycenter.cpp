// Geometric verification of barycenter placement: interfaces of symmetric
// inputs must lie on the bisecting plane, and 1-1-1-1 barycenters must land
// on the exact edge midpoints / face centers / tetrahedron center.

#include <delaunay_interfaces/interface_generation.hpp>
#include "test_helpers.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace delaunay_interfaces;
using test_helpers::approx_equal;
using test_helpers::count_simplices;

namespace {

bool point_approx_equal(const Point3D& p, double x, double y, double z) {
    return approx_equal(p.x(), x) && approx_equal(p.y(), y) && approx_equal(p.z(), z);
}

bool contains_point(const Points& points, double x, double y, double z) {
    for (const auto& p : points) {
        if (point_approx_equal(p, x, y, z)) return true;
    }
    return false;
}

// 3-1 coloring with the triangle at y=-1 and the apex at y=+1: every
// barycenter mixes equal numbers of y=-1 and y=+1 points, so the whole
// interface must lie on the y=0 plane.
void test_31_interface_lies_on_bisecting_plane() {
    std::cout << "Test: 3-1 interface lies on bisecting plane\n";

    double r = 1.0;
    Points points = {
        {r, -1.0, 0.0},
        {-r / 2.0, -1.0, r * std::sqrt(3.0) / 2.0},
        {-r / 2.0, -1.0, -r * std::sqrt(3.0) / 2.0},
        {0.0, 1.0, 0.0}
    };
    ColorLabels colors = {1, 1, 1, 2};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors);

    assert(!surface.vertices.empty());
    for (const auto& v : surface.vertices) {
        assert(approx_equal(v.y(), 0.0));
    }
    assert(count_simplices(surface.filtration.materialized()).triangles == 6);
}

// 2-2 coloring with color 1 at y=-1 and color 2 at y=+1: same symmetry
// argument, the interface must lie on y=0.
void test_22_interface_lies_on_bisecting_plane() {
    std::cout << "Test: 2-2 interface lies on bisecting plane\n";

    Points points = {
        {-1.0, -1.0, 0.0},
        {1.0, -1.0, 0.0},
        {0.0, 1.0, -1.0},
        {0.0, 1.0, 1.0}
    };
    ColorLabels colors = {1, 1, 2, 2};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors);

    assert(!surface.vertices.empty());
    for (const auto& v : surface.vertices) {
        assert(approx_equal(v.y(), 0.0));
    }
    assert(count_simplices(surface.filtration.materialized()).triangles == 8);
}

// Unit tetrahedron, all colors distinct: the 11 barycenters are exactly the
// 6 edge midpoints, 4 face centers, and the tetrahedron center.
void test_1111_barycenter_coordinates() {
    std::cout << "Test: 1-1-1-1 barycenter coordinates\n";

    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    ColorLabels colors = {1, 2, 3, 4};

    InterfaceGenerator generator;
    auto surface = generator.compute_interface_surface(points, colors);

    assert(surface.vertices.size() == 11);

    const double third = 1.0 / 3.0;
    const double expected[11][3] = {
        {0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5},   // edge midpoints
        {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
        {third, third, 0.0}, {third, 0.0, third},             // face centers
        {0.0, third, third}, {third, third, third},
        {0.25, 0.25, 0.25}                                    // tetrahedron center
    };
    for (const auto& e : expected) {
        assert(contains_point(surface.vertices, e[0], e[1], e[2]));
    }
}

} // namespace

int main() {
    test_31_interface_lies_on_bisecting_plane();
    test_22_interface_lies_on_bisecting_plane();
    test_1111_barycenter_coordinates();

    std::cout << "All barycenter verification tests passed\n";
    return 0;
}
