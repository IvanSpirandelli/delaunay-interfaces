#pragma once

#include "types.hpp"
#include <algorithm>
#include <map>

namespace delaunay_interfaces {

// Utility functions for chromatic partitioning
inline Partition get_chromatic_partitioning(
    const std::vector<int>& vertices,
    const ColorLabels& color_labels
) {
    std::map<int, std::vector<int>> parts_map;

    for (int vertex : vertices) {
        int color = color_labels[vertex];
        parts_map[color].push_back(vertex);
    }

    Partition parts;
    for (auto& [color, verts] : parts_map) {
        parts.push_back(verts);
    }

    // Sort by size (descending)
    std::sort(parts.begin(), parts.end(),
        [](const auto& a, const auto& b) { return a.size() > b.size(); }
    );

    return parts;
}

inline Partition get_chromatic_partitioning(
    const Tetrahedron& tet,
    const ColorLabels& color_labels
) {
    std::vector<int> vertices(tet.begin(), tet.end());
    return get_chromatic_partitioning(vertices, color_labels);
}

inline Point3D compute_barycenter(const Points& points, const std::vector<int>& indices) {
    Point3D center = Point3D::Zero();
    for (int idx : indices) {
        center += points[idx];
    }
    return center / static_cast<double>(indices.size());
}

inline double euclidean_distance(const Point3D& p1, const Point3D& p2) {
    return (p1 - p2).norm();
}

} // namespace delaunay_interfaces
