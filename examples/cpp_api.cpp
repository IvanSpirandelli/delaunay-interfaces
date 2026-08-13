#include <iostream>
#include <delaunay_interfaces/interface_generation.hpp>

using namespace delaunay_interfaces;

int main() {
    std::cout << "DelaunayInterfaces C++ - Simple Example\n";
    std::cout << "========================================\n\n";

    // Create a simple tetrahedral point cloud
    Points points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 1.0, 0.0},
        {0.5, 0.5, 1.0},
        {2.0, 0.0, 0.0},
        {2.5, 1.0, 0.0},
        {2.5, 0.5, 1.0},
        {1.5, 0.5, 0.5}
    };

    // Two groups of points (colors 1 and 2)
    ColorLabels colors = {1, 1, 1, 1, 2, 2, 2, 2};

    std::cout << "Input:\n";
    std::cout << "  Points: " << points.size() << "\n";
    std::cout << "  Colors: 2 (with " << colors.size() << " labels)\n\n";

    // Example 1: Regular Delaunay
    std::cout << "Example 1: Regular Delaunay Complex\n";
    std::cout << "------------------------------------\n";
    {
        InterfaceGenerator generator;
        auto surface = generator.compute_interface_surface(points, colors, {}, false, false);

        std::cout << "  Barycenters: " << surface.vertices.size() << "\n";
        std::cout << "  Filtration simplices: " << surface.filtration.size() << "\n\n";
    }

    // Example 2: Weighted Alpha Complex
    std::cout << "Example 2: Weighted Alpha Complex\n";
    std::cout << "----------------------------------\n";
    {
        Radii radii(points.size(), 0.3);

        InterfaceGenerator generator;
        auto surface = generator.compute_interface_surface(points, colors, radii, true, true);

        std::cout << "  Barycenters: " << surface.vertices.size() << "\n";
        std::cout << "  Filtration simplices: " << surface.filtration.size() << "\n\n";

        // Show some filtration values
        std::cout << "  First 10 filtration simplices:\n";
        for (size_t i = 0; i < std::min(size_t(10), surface.filtration.size()); ++i) {
            const auto& [simplex, value] = surface.filtration[i];
            std::cout << "    Simplex (dim=" << simplex.size()-1 << "): value=" << value << "\n";
        }
    }

    std::cout << "\nDone!\n";
    return 0;
}
