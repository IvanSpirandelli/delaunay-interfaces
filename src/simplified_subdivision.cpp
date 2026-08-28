#include "delaunay_interfaces/simplified_subdivision.hpp"
#include "subdivision_driver.hpp"
#include <algorithm>
#include <stdexcept>

namespace delaunay_interfaces {

SimplifiedSubdivision::SimplifiedSubdivision(
    const Points& points,
    const ColorLabels& color_labels
) : points_(points), color_labels_(color_labels) {}

int32_t SimplifiedSubdivision::get_or_create_vertex(int atom_a, int atom_b) {
    AtomKey key = {std::min(atom_a, atom_b), std::max(atom_a, atom_b)};

    auto it = vertex_map_.find(key);
    if (it != vertex_map_.end()) {
        return it->second.first;
    }

    int32_t id = next_vertex_id_++;
    Point3D midpoint = (points_[atom_a] + points_[atom_b]) / 2.0;
    double dist = (points_[atom_a] - points_[atom_b]).norm();

    midpoints_.push_back(midpoint);
    vertex_map_.emplace(key, std::make_pair(id, dist));
    return id;
}

void SimplifiedSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    process_simplex_impl(simplex_vertices.data(), static_cast<int>(simplex_vertices.size()));
}

void SimplifiedSubdivision::process_simplex_impl(const int* verts, int n_verts) {
    if (n_verts > 4) {
        throw std::invalid_argument(
            "process_simplex supports at most 4 vertices (edge, triangle, or tetrahedron)");
    }

    int part_colors[4];
    int part_atoms[4][4];
    int part_sizes[4];
    int n_parts = 0;
    for (int i = 0; i < n_verts; ++i) {
        const int v = verts[i];
        const int c = color_labels_[v];
        int p = 0;
        while (p < n_parts && part_colors[p] < c) ++p;
        if (p == n_parts || part_colors[p] != c) {
            for (int q = n_parts; q > p; --q) {
                part_colors[q] = part_colors[q - 1];
                part_sizes[q] = part_sizes[q - 1];
                std::copy(part_atoms[q - 1], part_atoms[q - 1] + 4, part_atoms[q]);
            }
            part_colors[p] = c;
            part_sizes[p] = 0;
            ++n_parts;
        }
        part_atoms[p][part_sizes[p]++] = v;
    }
    if (n_parts != 2) {
        throw std::invalid_argument(
            "SimplifiedSubdivision requires exactly 2 colors, got " +
            std::to_string(n_parts));
    }

    int a_part = 0, b_part = 1;
    if (part_sizes[1] > part_sizes[0]) std::swap(a_part, b_part);

    // At most a 2x2 cross product, so a fixed array suffices.
    int32_t cross_vertices[4];
    int n_cross = 0;
    for (int i = 0; i < part_sizes[a_part]; ++i) {
        for (int j = 0; j < part_sizes[b_part]; ++j) {
            cross_vertices[n_cross++] =
                get_or_create_vertex(part_atoms[a_part][i], part_atoms[b_part][j]);
        }
    }

    switch (n_cross) {
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
                               std::to_string(n_cross));
    }
}

void SimplifiedSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    process_simplex_impl(tet.data(), 4);
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
