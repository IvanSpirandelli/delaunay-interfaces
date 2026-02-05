#pragma once

#include <vector>
#include <tuple>
#include <cstdint>
#include <Eigen/Dense>

namespace delaunay_interfaces {

// Type aliases
using Point3D = Eigen::Vector3d;
using Points = std::vector<Point3D>;
using ColorLabels = std::vector<int>;
using Radii = std::vector<double>;
using Tetrahedron = std::array<int, 4>;
using Tetrahedra = std::vector<Tetrahedron>;
using Simplex = std::vector<int32_t>;
using SimplexWithFiltration = std::tuple<Simplex, double>;
using Filtration = std::vector<SimplexWithFiltration>;
using Partition = std::vector<std::vector<int>>;

// Configuration struct
struct ComplexConfig {
    bool weighted = true;
    bool alpha = true;

    ComplexConfig() = default;
    ComplexConfig(bool w, bool a) : weighted(w), alpha(a) {}
};

// Result structure
struct InterfaceSurface {
    Points vertices;
    Filtration filtration;
    bool weighted;
    bool alpha;
};

} // namespace delaunay_interfaces
