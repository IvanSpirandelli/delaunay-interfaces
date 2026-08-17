#pragma once

#include "types.hpp"
#include <algorithm>
#include <map>

namespace delaunay_interfaces {

// Groups vertex indices by color, parts sorted by descending size so the
// part order is deterministic.
[[nodiscard]] inline Partition compute_chromatic_partition(
    const std::vector<int>& vertices,
    const ColorLabels& color_labels
) {
    std::map<int, std::vector<int>> parts_map;
    for (int vertex : vertices) {
        parts_map[color_labels[vertex]].push_back(vertex);
    }

    Partition parts;
    parts.reserve(parts_map.size());
    for (auto& entry : parts_map) {
        parts.push_back(std::move(entry.second));
    }

    // Stable so that equal-size parts keep their color order (the map above
    // yields them sorted by color), making the output deterministic across
    // standard-library implementations.
    std::stable_sort(parts.begin(), parts.end(),
        [](const auto& a, const auto& b) { return a.size() > b.size(); }
    );

    return parts;
}

[[nodiscard]] inline Partition compute_chromatic_partition(
    const Tetrahedron& tet,
    const ColorLabels& color_labels
) {
    return compute_chromatic_partition({tet.begin(), tet.end()}, color_labels);
}

[[nodiscard]] inline Point3D compute_barycenter(const Points& points, const std::vector<int>& indices) {
    Point3D center = Point3D::Zero();
    for (int idx : indices) {
        center += points[idx];
    }
    return center / static_cast<double>(indices.size());
}

[[nodiscard]] inline double euclidean_distance(const Point3D& p1, const Point3D& p2) {
    return (p1 - p2).norm();
}

} // namespace delaunay_interfaces
