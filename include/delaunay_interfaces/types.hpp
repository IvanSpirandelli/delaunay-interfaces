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

// Simplex naming convention: simplices of the input complex (Delaunay /
// alpha) index input atoms as int and are named Tetrahedron / FreeTriangle /
// FreeEdge; simplices of the generated interface index the surface's own
// vertices as int32_t and are prefixed Surface. The interface within a
// generating d-simplex is (d-1)-dimensional, so it is a surface wherever the
// input complex is locally solid (sheets, joined along seams where >= 3
// colors meet) and degenerates with the complex where it thins out: free
// triangles yield curve pieces, free edges isolated vertices.

// Input complex.
using Tetrahedron = std::array<int, 4>;
using Tetrahedra = std::vector<Tetrahedron>;
using FreeTriangle = std::array<int, 3>;
using FreeEdge = std::array<int, 2>;
using Partition = std::vector<std::vector<int>>; // atom indices grouped by color

// Generated interface. Vertices are barycenters (barycentric construction)
// or midpoints (simplified construction). SurfaceSimplex has mixed
// dimension, hence a vector; the fixed-dimension faces are arrays.
using SurfaceSimplex = std::vector<int32_t>;
using FiltrationEntry = std::tuple<SurfaceSimplex, double>;
using Filtration = std::vector<FiltrationEntry>;
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
    MulticoloredSimplices generating_simplices;
    VertexAtomIndices vertex_atom_indices;
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
    MulticoloredSimplices generating_simplices; // Generating input simplices
    bool alpha;
};

} // namespace delaunay_interfaces
