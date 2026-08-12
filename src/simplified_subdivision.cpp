#include "delaunay_interfaces/simplified_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include <algorithm>
#include <stdexcept>

namespace delaunay_interfaces {

SimplifiedSubdivision::SimplifiedSubdivision(
    const Points& points,
    const ColorLabels& color_labels
) : points_(points), color_labels_(color_labels) {}

int32_t SimplifiedSubdivision::get_or_create_vertex(int atom_a, int atom_b) {
    std::array<int, 2> key = {std::min(atom_a, atom_b), std::max(atom_a, atom_b)};

    auto it = vertex_map_.find(key);
    if (it != vertex_map_.end()) {
        return it->second.first;
    }

    int32_t id = next_vertex_id_++;
    Point3D midpoint = (points_[atom_a] + points_[atom_b]) / 2.0;
    double dist = (points_[atom_a] - points_[atom_b]).norm();

    vertices_.push_back(midpoint);
    vertex_map_[key] = {id, dist};
    return id;
}

void SimplifiedSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    auto partition = get_chromatic_partitioning(simplex_vertices, color_labels_);
    if (partition.size() != 2) {
        throw std::invalid_argument(
            "SimplifiedSubdivision requires exactly 2 colors, got " +
            std::to_string(partition.size()));
    }

    const auto& part_a = partition[0]; // larger or equal group
    const auto& part_b = partition[1];

    // Enumerate all cross-color pairs → vertices
    std::vector<int32_t> cross_vertices;
    cross_vertices.reserve(part_a.size() * part_b.size());

    for (int a : part_a) {
        for (int b : part_b) {
            cross_vertices.push_back(get_or_create_vertex(a, b));
        }
    }

    size_t n = cross_vertices.size();

    if (n == 1) {
        // Free bicolored edge → single vertex (already registered)
        return;
    }

    if (n == 2) {
        // Free bicolored triangle (2-1 partition) → single edge
        std::array<int32_t, 2> edge = {cross_vertices[0], cross_vertices[1]};
        if (edge[0] > edge[1]) std::swap(edge[0], edge[1]);
        edges_.push_back(edge);
        return;
    }

    if (n == 3) {
        // 3-1 tetrahedron → single triangle
        Triangle tri = {cross_vertices[0], cross_vertices[1], cross_vertices[2]};
        std::sort(tri.begin(), tri.end());
        triangles_.push_back(tri);
        return;
    }

    if (n == 4) {
        // 2-2 tetrahedron → quad
        // cross_vertices layout (from nested loop): [a0b0, a0b1, a1b0, a1b1]
        // Cyclic order: a0b0 - a0b1 - a1b1 - a1b0
        //   (a0b0-a0b1 share a0, a0b1-a1b1 share b1, a1b1-a1b0 share a1, a1b0-a0b0 share b0)
        Quad quad = {cross_vertices[0], cross_vertices[1],
                     cross_vertices[3], cross_vertices[2]};
        quads_.push_back(quad);
        return;
    }

    // Should not reach here for a tetrahedron with 2 colors
    throw std::logic_error("Unexpected number of cross-color pairs: " + std::to_string(n));
}

void SimplifiedSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    std::vector<int> vertices(tet.begin(), tet.end());
    process_simplex(vertices);
}

std::vector<std::vector<int>> SimplifiedSubdivision::get_vertex_atom_indices() const {
    std::vector<std::vector<int>> result(next_vertex_id_);
    for (const auto& [key, id_val] : vertex_map_) {
        result[id_val.first] = {key[0], key[1]};
    }
    return result;
}

std::vector<double> SimplifiedSubdivision::get_vertex_filtration() const {
    std::vector<double> result(next_vertex_id_);
    for (const auto& [key, id_val] : vertex_map_) {
        result[id_val.first] = id_val.second;
    }
    return result;
}

} // namespace delaunay_interfaces
