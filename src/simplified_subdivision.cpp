#include "delaunay_interfaces/simplified_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "subdivision_driver.hpp"
#include <algorithm>
#include <stdexcept>

namespace delaunay_interfaces {

SimplifiedSubdivision::SimplifiedSubdivision(
    const Points& points,
    const ColorLabels& color_labels
) : points_(points), color_labels_(color_labels) {}

int32_t SimplifiedSubdivision::get_or_create_vertex(int atom_a, int atom_b) {
    std::array<int, 2> key = {std::min(atom_a, atom_b), std::max(atom_a, atom_b)};

    auto it = vertex_map_.lower_bound(key);
    if (it != vertex_map_.end() && it->first == key) {
        return it->second.first;
    }

    int32_t id = next_vertex_id_++;
    Point3D midpoint = (points_[atom_a] + points_[atom_b]) / 2.0;
    double dist = (points_[atom_a] - points_[atom_b]).norm();

    midpoints_.push_back(midpoint);
    vertex_map_.emplace_hint(it, key, std::make_pair(id, dist));
    return id;
}

void SimplifiedSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    auto partition = compute_chromatic_partition(simplex_vertices, color_labels_);
    if (partition.size() != 2) {
        throw std::invalid_argument(
            "SimplifiedSubdivision requires exactly 2 colors, got " +
            std::to_string(partition.size()));
    }

    // The part order from compute_chromatic_partition is deterministic; the
    // quad layout comes from the nested loop below and works for either order.
    const auto& part_a = partition[0];
    const auto& part_b = partition[1];

    std::vector<int32_t> cross_vertices;
    cross_vertices.reserve(part_a.size() * part_b.size());
    for (int a : part_a) {
        for (int b : part_b) {
            cross_vertices.push_back(get_or_create_vertex(a, b));
        }
    }

    switch (cross_vertices.size()) {
    case 1:
        // Free bicolored edge: single vertex, already registered above.
        return;
    case 2: {
        // Free bicolored triangle (2-1 partition): single edge.
        SurfaceEdge edge = {cross_vertices[0], cross_vertices[1]};
        if (edge[0] > edge[1]) std::swap(edge[0], edge[1]);
        edges_.push_back(edge);
        return;
    }
    case 3: {
        // 3-1 tetrahedron: single triangle.
        SurfaceTriangle tri = {cross_vertices[0], cross_vertices[1], cross_vertices[2]};
        std::sort(tri.begin(), tri.end());
        triangles_.push_back(tri);
        return;
    }
    case 4: {
        // 2-2 tetrahedron: quad. cross_vertices layout from the nested loop is
        // [a0b0, a0b1, a1b0, a1b1]; cyclic order a0b0 - a0b1 - a1b1 - a1b0
        // (consecutive quad vertices share an atom).
        SurfaceQuad quad = {cross_vertices[0], cross_vertices[1],
                     cross_vertices[3], cross_vertices[2]};
        quads_.push_back(quad);
        return;
    }
    default:
        throw std::logic_error("Unexpected number of cross-color pairs: " +
                               std::to_string(cross_vertices.size()));
    }
}

void SimplifiedSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    process_simplex({tet.begin(), tet.end()});
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

SimplifiedSurface compute_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha
) {
    detail::validate_inputs(points, color_labels, radii);

    const bool use_alpha = alpha.value_or(!radii.empty());
    SimplifiedSubdivision subdivision(points, color_labels);
    auto mc_simplices = detail::run_subdivision(subdivision, points, color_labels, radii, use_alpha);

    return SimplifiedSurface{
        subdivision.take_midpoints(),
        subdivision.take_triangles(),
        subdivision.take_quads(),
        subdivision.take_edges(),
        subdivision.get_vertex_atom_indices(),
        subdivision.get_vertex_filtration(),
        std::move(mc_simplices),
        use_alpha
    };
}

SimplifiedSurface compute_simplified_surface(
    const Points& points,
    const ColorLabels& color_labels,
    double radius
) {
    detail::validate_inputs(points, color_labels, Radii{});

    SimplifiedSubdivision subdivision(points, color_labels);
    auto mc_simplices = detail::run_subdivision(subdivision, points, color_labels, radius);

    return SimplifiedSurface{
        subdivision.take_midpoints(),
        subdivision.take_triangles(),
        subdivision.take_quads(),
        subdivision.take_edges(),
        subdivision.get_vertex_atom_indices(),
        subdivision.get_vertex_filtration(),
        std::move(mc_simplices),
        true
    };
}

} // namespace delaunay_interfaces
