#include "delaunay_interfaces/barycentric_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "delaunay_interfaces/interface_generation.hpp"
#include <algorithm>
#include <stdexcept>

namespace delaunay_interfaces {

BarycentricSubdivision::BarycentricSubdivision(
    const Points& points,
    const ColorLabels& color_labels,
    bool lower_star
) : points_(points), color_labels_(color_labels), lower_star_(lower_star) {}

double BarycentricSubdivision::star_value(double a, double b) const {
    return lower_star_ ? std::max(a, b) : std::min(a, b);
}

double BarycentricSubdivision::star_value(std::initializer_list<double> vals) const {
    return lower_star_ ? std::max(vals) : std::min(vals);
}

double BarycentricSubdivision::star_value(const std::vector<double>& vals) const {
    if (vals.empty()) return 0.0;
    double result = vals[0];
    for (size_t i = 1; i < vals.size(); ++i) {
        result = lower_star_ ? std::max(result, vals[i]) : std::min(result, vals[i]);
    }
    return result;
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
    size_t k = partitioning.size();
    if (k < 2) return 0.0;

    std::vector<Point3D> bcs;
    bcs.reserve(k);
    for (const auto& part : partitioning) {
        bcs.push_back(get_barycenter(part));
    }

    double sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            sum += euclidean_distance(bcs[i], bcs[j]);
            ++count;
        }
    }
    return sum / count;
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

void BarycentricSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    auto partition = get_chromatic_partitioning(simplex_vertices, color_labels_);
    size_t k = partition.size();
    if (k < 2) return;

    // --- 1. Enumerate all multicolored combinations ---
    // For each subset I ⊆ {0..k-1} with |I| >= 2,
    // enumerate all tuples of non-empty subsets S'_i ⊆ partition[i] for i ∈ I.

    struct Combination {
        std::vector<std::vector<int>> parts; // the partition for get_or_create_simplex
        std::vector<int> flat;               // sorted union of all vertices
    };

    std::vector<Combination> combinations;

    for (int mask = 3; mask < (1 << k); ++mask) {
        if (__builtin_popcount(mask) < 2) continue;

        // Collect part indices in this subset
        std::vector<int> part_indices;
        for (int i = 0; i < static_cast<int>(k); ++i) {
            if (mask & (1 << i)) part_indices.push_back(i);
        }

        // Cross-product of non-empty subsets of each selected part
        std::vector<std::vector<std::vector<int>>> partial = {{}};
        for (int pi : part_indices) {
            const auto& part = partition[pi];
            int n = static_cast<int>(part.size());
            std::vector<std::vector<std::vector<int>>> next;
            for (int smask = 1; smask < (1 << n); ++smask) {
                std::vector<int> subset;
                for (int j = 0; j < n; ++j) {
                    if (smask & (1 << j)) subset.push_back(part[j]);
                }
                for (const auto& prev : partial) {
                    auto extended = prev;
                    extended.push_back(subset);
                    next.push_back(std::move(extended));
                }
            }
            partial = std::move(next);
        }

        for (auto& tuple : partial) {
            Combination comb;
            comb.parts = std::move(tuple);
            for (const auto& s : comb.parts) {
                comb.flat.insert(comb.flat.end(), s.begin(), s.end());
            }
            std::sort(comb.flat.begin(), comb.flat.end());
            combinations.push_back(std::move(comb));
        }
    }

    size_t n = combinations.size();

    // --- 2. Get or create simplices, compute filtration values ---
    std::vector<std::pair<int32_t, double>> vertices(n);
    std::vector<bool> created(n);

    for (size_t i = 0; i < n; ++i) {
        auto info = get_or_create_simplex(combinations[i].parts);
        vertices[i] = {info.id, info.value};
        created[i] = info.newly_created;
    }

    // --- 3. Compute barycenters hierarchically ---
    // Min-level (smallest flat set size) → barycenter from original points
    // All others → average of min-level barycenters whose flat set ⊆ this flat set

    auto is_subset = [](const std::vector<int>& small, const std::vector<int>& big) {
        return std::includes(big.begin(), big.end(), small.begin(), small.end());
    };

    size_t min_level = combinations[0].flat.size();
    for (size_t i = 1; i < n; ++i) {
        min_level = std::min(min_level, combinations[i].flat.size());
    }

    std::vector<size_t> min_level_indices;
    for (size_t i = 0; i < n; ++i) {
        if (combinations[i].flat.size() == min_level) {
            min_level_indices.push_back(i);
        }
    }

    std::vector<Point3D> new_barycenters(n);

    // Min-level barycenters: flat average of original points
    for (size_t idx : min_level_indices) {
        new_barycenters[idx] = get_barycenter(combinations[idx].flat);
    }

    // All other barycenters: average of min-level barycenters with flat ⊆ this flat
    for (size_t i = 0; i < n; ++i) {
        if (combinations[i].flat.size() == min_level) continue;
        std::vector<Point3D> child_bcs;
        for (size_t midx : min_level_indices) {
            if (is_subset(combinations[midx].flat, combinations[i].flat)) {
                child_bcs.push_back(new_barycenters[midx]);
            }
        }
        new_barycenters[i] = get_barycenter_from_points(child_bcs);
    }

    // Add newly created barycenters
    for (size_t i = 0; i < n; ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // --- 4. Add edges (subset inclusion on flat sets) ---
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && is_subset(combinations[i].flat, combinations[j].flat)) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                double val = star_value(vertices[i].second, vertices[j].second);
                filtration_set_.insert({edge, val});
            }
        }
    }

    // --- 5. Add maximal chains as higher simplices ---
    // Group combinations by level (flat set size)
    std::map<size_t, std::vector<size_t>> levels;
    for (size_t i = 0; i < n; ++i) {
        levels[combinations[i].flat.size()].push_back(i);
    }

    std::vector<size_t> level_keys;
    for (const auto& [lv, _] : levels) {
        level_keys.push_back(lv);
    }

    if (level_keys.size() >= 3) {
        // Build chains: one element per level, each ⊂ next via flat inclusion
        std::vector<std::vector<size_t>> chains;
        for (size_t idx : levels[level_keys[0]]) {
            chains.push_back({idx});
        }

        for (size_t lv = 1; lv < level_keys.size(); ++lv) {
            std::vector<std::vector<size_t>> new_chains;
            for (const auto& chain : chains) {
                size_t last = chain.back();
                for (size_t idx : levels[level_keys[lv]]) {
                    if (is_subset(combinations[last].flat, combinations[idx].flat)) {
                        auto extended = chain;
                        extended.push_back(idx);
                        new_chains.push_back(std::move(extended));
                    }
                }
            }
            chains = std::move(new_chains);
        }

        // Each complete chain becomes a simplex (triangle for 3 levels, etc.)
        for (const auto& chain : chains) {
            Simplex simplex;
            std::vector<double> vals;
            for (size_t idx : chain) {
                simplex.push_back(vertices[idx].first);
                vals.push_back(vertices[idx].second);
            }
            std::sort(simplex.begin(), simplex.end());
            filtration_set_.insert({simplex, star_value(vals)});
        }
    }

    // --- 6. Add vertices to filtration ---
    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

void BarycentricSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    std::vector<int> vertices(tet.begin(), tet.end());
    process_simplex(vertices);
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

std::tuple<Points, Filtration, MulticoloredSimplices> get_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha,
    bool lower_star
) {
    if (points.size() != color_labels.size()) {
        throw std::invalid_argument("Each point must have a corresponding color_label");
    }

    if (weighted && radii.size() != points.size()) {
        throw std::invalid_argument("Each point must have an assigned radius for weighted complexes");
    }

    InterfaceGenerator generator;
    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    MulticoloredSimplices mc_simplices;

    if (weighted && alpha) {
        mc_simplices = generator.get_multicolored_simplices_weighted_alpha(
            points, color_labels, radii);

        for (const auto& tet : mc_simplices.generating_tetrahedra) {
            subdivision.process_tetrahedron(tet);
        }
        for (const auto& tri : mc_simplices.generating_free_triangles) {
            subdivision.process_simplex(tri);
        }
        for (const auto& edge : mc_simplices.generating_free_edges) {
            subdivision.process_simplex(edge);
        }
    } else {
        auto tetrahedra = generator.get_multicolored_tetrahedra(
            points, color_labels, radii, weighted, alpha);

        mc_simplices.generating_tetrahedra = tetrahedra;

        for (const auto& tet : tetrahedra) {
            subdivision.process_tetrahedron(tet);
        }
    }

    return {subdivision.get_barycenters(), subdivision.get_filtration(), std::move(mc_simplices)};
}

} // namespace delaunay_interfaces
