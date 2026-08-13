#pragma once

#include <array>
#include <vector>
#include <tuple>
#include <cstdint>
#include <Eigen/Dense>

namespace delaunay_interfaces {

using Point3D = Eigen::Vector3d;
using Points = std::vector<Point3D>;
using ColorLabels = std::vector<int>;
using Radii = std::vector<double>;
using Tetrahedron = std::array<int, 4>;
using Tetrahedra = std::vector<Tetrahedron>;
using FreeTriangle = std::array<int, 3>;
using FreeEdge = std::array<int, 2>;
using Simplex = std::vector<int32_t>;
using SimplexWithFiltration = std::tuple<Simplex, double>;
using Filtration = std::vector<SimplexWithFiltration>;
using Partition = std::vector<std::vector<int>>;

// Multicolored simplices from the complex (generating tetrahedra + free lower-dimensional)
struct MulticoloredSimplices {
    Tetrahedra generating_tetrahedra;
    std::vector<FreeTriangle> generating_free_triangles;
    std::vector<FreeEdge> generating_free_edges;
};

// Vertex-to-atom mapping: vertex_atom_indices[i] = sorted input atom indices for vertex i
using VertexAtomIndices = std::vector<std::vector<int>>;

struct InterfaceSurface {
    Points vertices;
    Filtration filtration;
    MulticoloredSimplices simplices;
    VertexAtomIndices vertex_atom_indices;
    bool weighted;
    bool alpha;
    bool lower_star;
};

// Simplified 2-color surface: each vertex = one cross-color atom pair midpoint
using Triangle = std::array<int32_t, 3>;
using Quad = std::array<int32_t, 4>;

struct SimplifiedSurface {
    Points vertices;                         // Midpoints of cross-color edges
    std::vector<Triangle> triangles;         // Triangle faces (from 3-1 tets)
    std::vector<Quad> quads;                 // Quad faces (from 2-2 tets)
    std::vector<std::array<int32_t, 2>> edges;  // Free edges (from free triangles)
    VertexAtomIndices vertex_atom_indices;    // Each entry has exactly 2 atom indices
    std::vector<double> vertex_filtration;   // Filtration value per vertex
    MulticoloredSimplices simplices;         // Generating structures
    bool weighted;
    bool alpha;
};

} // namespace delaunay_interfaces
