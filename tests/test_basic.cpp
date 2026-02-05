#include <iostream>
#include <cassert>
#include <cmath>
#include <delaunay_interfaces/interface_generation.hpp>
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

int main() {
    std::cout << "Running DelaunayInterfaces C++ Tests\n";
    std::cout << "=====================================\n\n";

    try {
        test_chromatic_partitioning();
        test_barycenter();
        test_simple_delaunay();
        test_weighted_alpha();
        test_input_validation();

        std::cout << "\nAll tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << "\n";
        return 1;
    }
}
