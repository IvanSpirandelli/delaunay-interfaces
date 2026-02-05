#include "delaunay_interfaces/barycentric_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "delaunay_interfaces/interface_generation.hpp"
#include <algorithm>
#include <stdexcept>

namespace delaunay_interfaces {

BarycentricSubdivision::BarycentricSubdivision(
    const Points& points,
    const ColorLabels& color_labels
) : points_(points), color_labels_(color_labels) {}

Partition BarycentricSubdivision::get_chromatic_partitioning(const Tetrahedron& tet) const {
    return delaunay_interfaces::get_chromatic_partitioning(tet, color_labels_);
}

Point3D BarycentricSubdivision::get_barycenter(const std::vector<int>& vertices) const {
    return compute_barycenter(points_, vertices);
}

Point3D BarycentricSubdivision::get_barycenter_from_points(const std::vector<Point3D>& points) const {
    Point3D center = Point3D::Zero();
    for (const auto& p : points) {
        center += p;
    }
    return center / static_cast<double>(points.size());
}

double BarycentricSubdivision::compute_filtration_value(const Partition& partitioning) const {
    if (partitioning.size() == 2) {
        auto bc1 = get_barycenter(partitioning[0]);
        auto bc2 = get_barycenter(partitioning[1]);
        return euclidean_distance(bc1, bc2);
    } else if (partitioning.size() == 3) {
        auto bc1 = get_barycenter(partitioning[0]);
        auto bc2 = get_barycenter(partitioning[1]);
        auto bc3 = get_barycenter(partitioning[2]);
        double a = euclidean_distance(bc1, bc2);
        double b = euclidean_distance(bc1, bc3);
        double c = euclidean_distance(bc2, bc3);
        return (a + b + c) / 3.0;
    } else if (partitioning.size() == 4) {
        auto bc1 = get_barycenter(partitioning[0]);
        auto bc2 = get_barycenter(partitioning[1]);
        auto bc3 = get_barycenter(partitioning[2]);
        auto bc4 = get_barycenter(partitioning[3]);
        double a = euclidean_distance(bc1, bc2);
        double b = euclidean_distance(bc1, bc3);
        double c = euclidean_distance(bc1, bc4);
        double d = euclidean_distance(bc2, bc3);
        double e = euclidean_distance(bc2, bc4);
        double f = euclidean_distance(bc3, bc4);
        return (a + b + c + d + e + f) / 6.0;
    }
    return 0.0;
}

BarycentricSubdivision::SimplexInfo BarycentricSubdivision::get_or_create_simplex(
    const std::vector<std::vector<int>>& partitioning
) {
    // Create sorted key from all vertices
    std::vector<int> key;
    for (const auto& part : partitioning) {
        key.insert(key.end(), part.begin(), part.end());
    }
    std::sort(key.begin(), key.end());

    auto it = simplex_map_.find(key);
    if (it != simplex_map_.end()) {
        return SimplexInfo{it->second.first, it->second.second, false};
    } else {
        int32_t id = next_simplex_id_++;
        double value = compute_filtration_value(partitioning);
        simplex_map_[key] = {id, value};
        return SimplexInfo{id, value, true};
    }
}

void BarycentricSubdivision::extend_scaffold_2_2(
    const std::vector<int>& part1,
    const std::vector<int>& part2
) {
    // 2-2 partitioning: [u,v] vs [x,y]
    int u = part1[0], v = part1[1];
    int x = part2[0], y = part2[1];

    // Define the 9 multicolored combinations
    // Indices 0-3: edges (2 vertices), 4-7: triangles (3 vertices), 8: tet (4 vertices)
    std::vector<std::vector<std::vector<int>>> mc_combinations = {
        {{u},{x}}, {{v},{x}}, {{v},{y}}, {{u},{y}},      // 0-3: edges
        {{u,v},{x}}, {{v},{x,y}}, {{u,v},{y}}, {{u},{x,y}},  // 4-7: triangles
        {{u,v},{x,y}}  // 8: tet
    };

    std::vector<std::pair<int32_t, double>> vertices;
    std::vector<bool> created;

    // Create or get simplex IDs
    for (const auto& comb : mc_combinations) {
        auto info = get_or_create_simplex(comb);
        vertices.push_back({info.id, info.value});
        created.push_back(info.newly_created);
    }

    // Compute barycenters hierarchically (matching Julia implementation)
    std::vector<Point3D> new_barycenters(9);

    // Step 1: Compute edge barycenters (indices 0-3) from original points
    std::vector<int> edge_indices_list = {0, 1, 2, 3};
    for (int idx : edge_indices_list) {
        std::vector<int> all_verts;
        for (const auto& part : mc_combinations[idx]) {
            all_verts.insert(all_verts.end(), part.begin(), part.end());
        }
        new_barycenters[idx] = get_barycenter(all_verts);
    }

    // Helper to flatten mc_combination to vertex set
    auto flatten = [](const std::vector<std::vector<int>>& comb) {
        std::vector<int> result;
        for (const auto& part : comb) {
            result.insert(result.end(), part.begin(), part.end());
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Helper to check if one set is subset of another
    auto is_subset = [](const std::vector<int>& small, const std::vector<int>& big) {
        return std::includes(big.begin(), big.end(), small.begin(), small.end());
    };

    // Flatten all combinations for subset checks
    std::vector<std::vector<int>> flat_combs(9);
    for (size_t i = 0; i < mc_combinations.size(); ++i) {
        flat_combs[i] = flatten(mc_combinations[i]);
    }

    // Step 2: Compute triangle barycenters (indices 4-7) from edge barycenters
    std::vector<int> tri_indices_list = {4, 5, 6, 7};
    for (int tri_idx : tri_indices_list) {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[tri_idx])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[tri_idx] = get_barycenter_from_points(edge_bcs);
    }

    // Step 3: Compute tet barycenter (index 8) from edge barycenters
    {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[8])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[8] = get_barycenter_from_points(edge_bcs);
    }

    // Add newly created barycenters
    for (size_t i = 0; i < created.size(); ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // Add edges dynamically based on subset relationships (matching Julia implementation)
    for (size_t i = 0; i < flat_combs.size(); ++i) {
        for (size_t j = 0; j < flat_combs.size(); ++j) {
            if (i != j && is_subset(flat_combs[i], flat_combs[j])) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                double val = std::min(vertices[i].second, vertices[j].second);
                filtration_set_.insert({edge, val});
            }
        }
    }

    // Add triangles (8 triangles as in Julia code)
    std::vector<std::tuple<int, int, int>> triangle_indices = {
        {8, 0, 4}, {8, 4, 1}, {8, 1, 5}, {8, 5, 2},
        {8, 2, 6}, {8, 6, 3}, {8, 3, 7}, {8, 7, 0}
    };

    for (const auto& [i, j, k] : triangle_indices) {
        Simplex tri = {vertices[i].first, vertices[j].first, vertices[k].first};
        std::sort(tri.begin(), tri.end());
        double val = std::min({vertices[i].second, vertices[j].second, vertices[k].second});
        filtration_set_.insert({tri, val});
    }

    // Add vertices to filtration
    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

void BarycentricSubdivision::extend_scaffold_3_1(
    const std::vector<int>& part1,
    const std::vector<int>& part2
) {
    // 3-1 partitioning: [u,v,w] vs [x]
    int u = part1[0], v = part1[1], w = part1[2];
    int x = part2[0];

    // Indices 0-2: edges (2 vertices), 3-5: triangles (3 vertices), 6: tet (4 vertices)
    std::vector<std::vector<std::vector<int>>> mc_combinations = {
        {{u},{x}}, {{v},{x}}, {{w},{x}},        // 0-2: edges
        {{u,v},{x}}, {{v,w},{x}}, {{u,w},{x}},  // 3-5: triangles
        {{u,v,w},{x}}  // 6: tet
    };

    std::vector<std::pair<int32_t, double>> vertices;
    std::vector<bool> created;

    for (const auto& comb : mc_combinations) {
        auto info = get_or_create_simplex(comb);
        vertices.push_back({info.id, info.value});
        created.push_back(info.newly_created);
    }

    // Compute barycenters hierarchically
    std::vector<Point3D> new_barycenters(7);

    // Helper to flatten mc_combination to vertex set
    auto flatten = [](const std::vector<std::vector<int>>& comb) {
        std::vector<int> result;
        for (const auto& part : comb) {
            result.insert(result.end(), part.begin(), part.end());
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Helper to check if one set is subset of another
    auto is_subset = [](const std::vector<int>& small, const std::vector<int>& big) {
        return std::includes(big.begin(), big.end(), small.begin(), small.end());
    };

    // Flatten all combinations
    std::vector<std::vector<int>> flat_combs(7);
    for (size_t i = 0; i < mc_combinations.size(); ++i) {
        flat_combs[i] = flatten(mc_combinations[i]);
    }

    // Step 1: Compute edge barycenters (indices 0-2) from original points
    std::vector<int> edge_indices_list = {0, 1, 2};
    for (int idx : edge_indices_list) {
        std::vector<int> all_verts;
        for (const auto& part : mc_combinations[idx]) {
            all_verts.insert(all_verts.end(), part.begin(), part.end());
        }
        new_barycenters[idx] = get_barycenter(all_verts);
    }

    // Step 2: Compute triangle barycenters (indices 3-5) from edge barycenters
    std::vector<int> tri_indices_list = {3, 4, 5};
    for (int tri_idx : tri_indices_list) {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[tri_idx])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[tri_idx] = get_barycenter_from_points(edge_bcs);
    }

    // Step 3: Compute tet barycenter (index 6) from edge barycenters
    {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[6])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[6] = get_barycenter_from_points(edge_bcs);
    }

    for (size_t i = 0; i < created.size(); ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // Add edges dynamically based on subset relationships (matching Julia implementation)
    for (size_t i = 0; i < flat_combs.size(); ++i) {
        for (size_t j = 0; j < flat_combs.size(); ++j) {
            if (i != j && is_subset(flat_combs[i], flat_combs[j])) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                double val = std::min(vertices[i].second, vertices[j].second);
                filtration_set_.insert({edge, val});
            }
        }
    }

    // Add triangles (6 triangles)
    std::vector<std::tuple<int, int, int>> triangle_indices = {
        {6, 0, 3}, {6, 3, 1}, {6, 1, 4},
        {6, 4, 2}, {6, 2, 5}, {6, 5, 0}
    };

    for (const auto& [i, j, k] : triangle_indices) {
        Simplex tri = {vertices[i].first, vertices[j].first, vertices[k].first};
        std::sort(tri.begin(), tri.end());
        double val = std::min({vertices[i].second, vertices[j].second, vertices[k].second});
        filtration_set_.insert({tri, val});
    }

    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

void BarycentricSubdivision::extend_scaffold_2_1_1(
    const std::vector<int>& part1,
    const std::vector<int>& part2,
    const std::vector<int>& part3
) {
    // 2-1-1 partitioning: [a,b] vs [u] vs [x]
    int a = part1[0], b = part1[1];
    int u = part2[0];
    int x = part3[0];

    // Indices 0-4: edges (2 vertices), 5-8: triangles (3 vertices), 9: tet (4 vertices)
    std::vector<std::vector<std::vector<int>>> mc_combinations = {
        {{a},{u}}, {{a},{x}}, {{b},{x}}, {{b},{u}},  // 0-3: edges
        {{u},{x}},  // 4: edge
        {{a},{u},{x}}, {{a,b},{x}}, {{b},{u},{x}}, {{a,b},{u}},  // 5-8: triangles
        {{a,b},{u},{x}}  // 9: tet
    };

    std::vector<std::pair<int32_t, double>> vertices;
    std::vector<bool> created;

    for (const auto& comb : mc_combinations) {
        auto info = get_or_create_simplex(comb);
        vertices.push_back({info.id, info.value});
        created.push_back(info.newly_created);
    }

    // Compute barycenters hierarchically
    std::vector<Point3D> new_barycenters(10);

    // Helper to flatten mc_combination to vertex set
    auto flatten = [](const std::vector<std::vector<int>>& comb) {
        std::vector<int> result;
        for (const auto& part : comb) {
            result.insert(result.end(), part.begin(), part.end());
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Helper to check if one set is subset of another
    auto is_subset = [](const std::vector<int>& small, const std::vector<int>& big) {
        return std::includes(big.begin(), big.end(), small.begin(), small.end());
    };

    // Flatten all combinations
    std::vector<std::vector<int>> flat_combs(10);
    for (size_t i = 0; i < mc_combinations.size(); ++i) {
        flat_combs[i] = flatten(mc_combinations[i]);
    }

    // Step 1: Compute edge barycenters (indices 0-4) from original points
    std::vector<int> edge_indices_list = {0, 1, 2, 3, 4};
    for (int idx : edge_indices_list) {
        std::vector<int> all_verts;
        for (const auto& part : mc_combinations[idx]) {
            all_verts.insert(all_verts.end(), part.begin(), part.end());
        }
        new_barycenters[idx] = get_barycenter(all_verts);
    }

    // Step 2: Compute triangle barycenters (indices 5-8) from edge barycenters
    std::vector<int> tri_indices_list = {5, 6, 7, 8};
    for (int tri_idx : tri_indices_list) {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[tri_idx])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[tri_idx] = get_barycenter_from_points(edge_bcs);
    }

    // Step 3: Compute tet barycenter (index 9) from edge barycenters
    {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[9])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[9] = get_barycenter_from_points(edge_bcs);
    }

    for (size_t i = 0; i < created.size(); ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // Add edges dynamically based on subset relationships (matching Julia implementation)
    for (size_t i = 0; i < flat_combs.size(); ++i) {
        for (size_t j = 0; j < flat_combs.size(); ++j) {
            if (i != j && is_subset(flat_combs[i], flat_combs[j])) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                double val = std::min(vertices[i].second, vertices[j].second);
                filtration_set_.insert({edge, val});
            }
        }
    }

    // Add triangles (10 triangles)
    std::vector<std::tuple<int, int, int>> triangle_indices = {
        {9, 0, 5}, {9, 5, 4}, {9, 4, 7}, {9, 7, 3},
        {9, 3, 8}, {9, 8, 0}, {9, 2, 7}, {9, 5, 1},
        {9, 1, 6}, {9, 6, 2}
    };

    for (const auto& [i, j, k] : triangle_indices) {
        Simplex tri = {vertices[i].first, vertices[j].first, vertices[k].first};
        std::sort(tri.begin(), tri.end());
        double val = std::min({vertices[i].second, vertices[j].second, vertices[k].second});
        filtration_set_.insert({tri, val});
    }

    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

void BarycentricSubdivision::extend_scaffold_1_1_1_1(
    const std::vector<int>& part1,
    const std::vector<int>& part2,
    const std::vector<int>& part3,
    const std::vector<int>& part4
) {
    // 1-1-1-1 partitioning: [a] vs [i] vs [u] vs [x]
    int a = part1[0];
    int i_idx = part2[0];
    int u = part3[0];
    int x = part4[0];

    // Indices 0-5: edges (2 vertices), 6-9: triangles (3 vertices), 10: tet (4 vertices)
    std::vector<std::vector<std::vector<int>>> mc_combinations = {
        {{a},{i_idx}}, {{a},{u}}, {{a},{x}},  // 0-2: edges
        {{i_idx},{u}}, {{i_idx},{x}}, {{u},{x}},  // 3-5: edges
        {{a},{i_idx},{u}}, {{a},{i_idx},{x}}, {{i_idx},{u},{x}}, {{a},{u},{x}},  // 6-9: triangles
        {{a},{i_idx},{u},{x}}  // 10: tet
    };

    std::vector<std::pair<int32_t, double>> vertices;
    std::vector<bool> created;

    for (const auto& comb : mc_combinations) {
        auto info = get_or_create_simplex(comb);
        vertices.push_back({info.id, info.value});
        created.push_back(info.newly_created);
    }

    // Compute barycenters hierarchically
    std::vector<Point3D> new_barycenters(11);

    // Helper to flatten mc_combination to vertex set
    auto flatten = [](const std::vector<std::vector<int>>& comb) {
        std::vector<int> result;
        for (const auto& part : comb) {
            result.insert(result.end(), part.begin(), part.end());
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Helper to check if one set is subset of another
    auto is_subset = [](const std::vector<int>& small, const std::vector<int>& big) {
        return std::includes(big.begin(), big.end(), small.begin(), small.end());
    };

    // Flatten all combinations
    std::vector<std::vector<int>> flat_combs(11);
    for (size_t i = 0; i < mc_combinations.size(); ++i) {
        flat_combs[i] = flatten(mc_combinations[i]);
    }

    // Step 1: Compute edge barycenters (indices 0-5) from original points
    std::vector<int> edge_indices_list = {0, 1, 2, 3, 4, 5};
    for (int idx : edge_indices_list) {
        std::vector<int> all_verts;
        for (const auto& part : mc_combinations[idx]) {
            all_verts.insert(all_verts.end(), part.begin(), part.end());
        }
        new_barycenters[idx] = get_barycenter(all_verts);
    }

    // Step 2: Compute triangle barycenters (indices 6-9) from edge barycenters
    std::vector<int> tri_indices_list = {6, 7, 8, 9};
    for (int tri_idx : tri_indices_list) {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[tri_idx])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[tri_idx] = get_barycenter_from_points(edge_bcs);
    }

    // Step 3: Compute tet barycenter (index 10) from edge barycenters
    {
        std::vector<Point3D> edge_bcs;
        for (int edge_idx : edge_indices_list) {
            if (is_subset(flat_combs[edge_idx], flat_combs[10])) {
                edge_bcs.push_back(new_barycenters[edge_idx]);
            }
        }
        new_barycenters[10] = get_barycenter_from_points(edge_bcs);
    }

    for (size_t i = 0; i < created.size(); ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // Add edges dynamically based on subset relationships (matching Julia implementation)
    for (size_t i = 0; i < flat_combs.size(); ++i) {
        for (size_t j = 0; j < flat_combs.size(); ++j) {
            if (i != j && is_subset(flat_combs[i], flat_combs[j])) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                double val = std::min(vertices[i].second, vertices[j].second);
                filtration_set_.insert({edge, val});
            }
        }
    }

    // Add triangles (12 triangles)
    std::vector<std::tuple<int, int, int>> triangle_indices = {
        {10, 3, 8}, {10, 8, 4}, {10, 4, 7}, {10, 7, 0},
        {10, 0, 6}, {10, 6, 3}, {10, 9, 1}, {10, 5, 9},
        {10, 8, 5}, {10, 1, 6}, {10, 9, 2}, {10, 2, 7}
    };

    for (const auto& [i, j, k] : triangle_indices) {
        Simplex tri = {vertices[i].first, vertices[j].first, vertices[k].first};
        std::sort(tri.begin(), tri.end());
        double val = std::min({vertices[i].second, vertices[j].second, vertices[k].second});
        filtration_set_.insert({tri, val});
    }

    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

void BarycentricSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    auto parts = get_chromatic_partitioning(tet);

    if (parts.size() == 2) {
        if (parts[0].size() == 2 && parts[1].size() == 2) {
            extend_scaffold_2_2(parts[0], parts[1]);
        } else if (parts[0].size() == 3 && parts[1].size() == 1) {
            extend_scaffold_3_1(parts[0], parts[1]);
        } else if (parts[0].size() == 1 && parts[1].size() == 3) {
            extend_scaffold_3_1(parts[1], parts[0]);
        } else {
            throw std::runtime_error("Invalid 2-part partitioning");
        }
    } else if (parts.size() == 3) {
        extend_scaffold_2_1_1(parts[0], parts[1], parts[2]);
    } else if (parts.size() == 4) {
        extend_scaffold_1_1_1_1(parts[0], parts[1], parts[2], parts[3]);
    }
}

Filtration BarycentricSubdivision::get_filtration() const {
    Filtration result(filtration_set_.begin(), filtration_set_.end());

    // Sort by simplex size, then by filtration value
    std::sort(result.begin(), result.end(),
        [](const auto& a, const auto& b) {
            if (std::get<0>(a).size() != std::get<0>(b).size()) {
                return std::get<0>(a).size() < std::get<0>(b).size();
            }
            return std::get<1>(a) < std::get<1>(b);
        }
    );

    return result;
}

std::pair<Points, Filtration> get_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha
) {
    if (points.size() != color_labels.size()) {
        throw std::invalid_argument("Each point must have a corresponding color_label");
    }

    if (weighted && radii.size() != points.size()) {
        throw std::invalid_argument("Each point must have an assigned radius for weighted complexes");
    }

    InterfaceGenerator generator;
    auto tetrahedra = generator.get_multicolored_tetrahedra(points, color_labels, radii, weighted, alpha);

    BarycentricSubdivision subdivision(points, color_labels);

    for (const auto& tet : tetrahedra) {
        subdivision.process_tetrahedron(tet);
    }

    return {subdivision.get_barycenters(), subdivision.get_filtration()};
}

} // namespace delaunay_interfaces
