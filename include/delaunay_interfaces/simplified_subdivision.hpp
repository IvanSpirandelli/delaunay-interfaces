#pragma once

#include "types.hpp"
#include <array>
#include <cstdint>
#include <unordered_map>

namespace delaunay_interfaces {

// Simplified subdivision for 2-color interfaces.
// Each vertex = midpoint of a cross-color edge (exactly 2 atom indices).
// Faces: 3-1 tet → triangle, 2-2 tet → quad.
class SimplifiedSubdivision {
public:
    // points and color_labels are held by reference and must outlive this object.
    SimplifiedSubdivision(const Points& points, const ColorLabels& color_labels);

    void process_tetrahedron(const Tetrahedron& tet);
    void process_simplex(const std::vector<int>& simplex_vertices);

    // Capacity hint from the shared subdivision driver; the counts it passes
    // are the barycentric construction's, so this is a no-op here.
    void reserve(size_t /*vertices*/, size_t /*edges*/, size_t /*chains*/) {}

    [[nodiscard]] const Points& get_midpoints() const { return midpoints_; }
    [[nodiscard]] const std::vector<SurfaceTriangle>& get_triangles() const { return triangles_; }
    [[nodiscard]] const std::vector<SurfaceQuad>& get_quads() const { return quads_; }
    [[nodiscard]] const std::vector<SurfaceEdge>& get_edges() const { return edges_; }
    [[nodiscard]] std::vector<std::vector<int>> get_vertex_atom_indices() const;
    [[nodiscard]] std::vector<double> get_vertex_filtration() const;

    // Move-out variants for result assembly; the moved-from field is spent.
    [[nodiscard]] Points take_midpoints() { return std::move(midpoints_); }
    [[nodiscard]] std::vector<SurfaceTriangle> take_triangles() { return std::move(triangles_); }
    [[nodiscard]] std::vector<SurfaceQuad> take_quads() { return std::move(quads_); }
    [[nodiscard]] std::vector<SurfaceEdge> take_edges() { return std::move(edges_); }

private:
    // Shared allocation-free core of process_simplex / process_tetrahedron.
    void process_simplex_impl(const int* verts, int n_verts);

    int32_t get_or_create_vertex(int atom_a, int atom_b);

    const Points& points_;
    const ColorLabels& color_labels_;

    Points midpoints_;
    std::vector<SurfaceTriangle> triangles_;
    std::vector<SurfaceQuad> quads_;
    std::vector<SurfaceEdge> edges_;

    // Sorted atom pair packed into a fixed-size key. Every vertex is the
    // midpoint of one cross-color edge, so an atom pair always fits.
    using AtomKey = std::array<int, 2>;
    struct AtomKeyHash {
        size_t operator()(const AtomKey& k) const {
            uint64_t h = 0;
            for (int a : k) {
                // splitmix64 round per element
                uint64_t x = h ^ (static_cast<uint64_t>(static_cast<uint32_t>(a)) + 0x9e3779b97f4a7c15ULL);
                x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
                x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
                h = x ^ (x >> 31);
            }
            return static_cast<size_t>(h);
        }
    };

    // Map from packed sorted atom pairs to (vertex_id, filtration_value)
    std::unordered_map<AtomKey, std::pair<int32_t, double>, AtomKeyHash> vertex_map_;
    int32_t next_vertex_id_ = 0;
};

} // namespace delaunay_interfaces
