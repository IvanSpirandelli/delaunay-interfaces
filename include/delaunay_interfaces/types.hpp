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

// Simplices of the input complex (Delaunay / alpha), indexing input atoms.
using Tetrahedron = std::array<int, 4>;
using Tetrahedra = std::vector<Tetrahedron>;
using FreeTriangle = std::array<int, 3>;
using FreeEdge = std::array<int, 2>;
using Partition = std::vector<std::vector<int>>; // atom indices grouped by color

// Simplices of the generated interface, indexing its own (barycenter or
// midpoint) vertices. Filtration simplices have mixed dimension, hence a
// vector; the simplified-surface faces have fixed dimension.
using Simplex = std::vector<int32_t>;
using SimplexWithFiltration = std::tuple<Simplex, double>;
using Filtration = std::vector<SimplexWithFiltration>;
using SurfaceTriangle = std::array<int32_t, 3>;
using SurfaceQuad = std::array<int32_t, 4>;
using SurfaceEdge = std::array<int32_t, 2>;

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
struct SimplifiedSurface {
    Points vertices;                         // Midpoints of cross-color edges
    std::vector<SurfaceTriangle> triangles;  // Triangle faces (from 3-1 tets)
    std::vector<SurfaceQuad> quads;          // Quad faces (from 2-2 tets)
    std::vector<SurfaceEdge> edges;          // Edges (one per free 2-1 input triangle)
    VertexAtomIndices vertex_atom_indices;    // Each entry has exactly 2 atom indices
    std::vector<double> vertex_filtration;   // Filtration value per vertex
    MulticoloredSimplices simplices;         // Generating structures
    bool weighted;
    bool alpha;
};

} // namespace delaunay_interfaces
