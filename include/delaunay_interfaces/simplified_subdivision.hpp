#pragma once

#include "types.hpp"
#include <map>

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

    [[nodiscard]] const Points& get_vertices() const { return vertices_; }
    [[nodiscard]] const std::vector<Triangle>& get_triangles() const { return triangles_; }
    [[nodiscard]] const std::vector<Quad>& get_quads() const { return quads_; }
    [[nodiscard]] const std::vector<std::array<int32_t, 2>>& get_edges() const { return edges_; }
    [[nodiscard]] std::vector<std::vector<int>> get_vertex_atom_indices() const;
    [[nodiscard]] std::vector<double> get_vertex_filtration() const;

private:
    int32_t get_or_create_vertex(int atom_a, int atom_b);

    const Points& points_;
    const ColorLabels& color_labels_;

    Points vertices_;
    std::vector<Triangle> triangles_;
    std::vector<Quad> quads_;
    std::vector<std::array<int32_t, 2>> edges_;

    // Map from sorted atom pair → (vertex_id, filtration_value)
    std::map<std::array<int, 2>, std::pair<int32_t, double>> vertex_map_;
    int32_t next_vertex_id_ = 0;
};

} // namespace delaunay_interfaces
