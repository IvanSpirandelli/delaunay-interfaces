#include "delaunay_interfaces/barycentric_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "subdivision_driver.hpp"
#include <algorithm>
#include <bitset>

namespace delaunay_interfaces {

namespace {

// One multicolored sub-simplex of the processed simplex: a choice of >= 2
// parts and a non-empty subset of each chosen part.
struct Combination {
    std::vector<std::vector<int>> parts; // the sub-partition, one entry per chosen part
    std::vector<int> flat;               // sorted union of all chosen vertices
};

// Enumerates every subset I of the parts with |I| >= 2, crossed with every
// tuple of non-empty subsets of the parts in I.
std::vector<Combination> enumerate_multicolored_combinations(const Partition& partition) {
    const int k = static_cast<int>(partition.size());
    std::vector<Combination> combinations;

    for (int mask = 0; mask < (1 << k); ++mask) {
        if (std::bitset<32>(mask).count() < 2) continue;

        std::vector<int> part_indices;
        for (int i = 0; i < k; ++i) {
            if (mask & (1 << i)) part_indices.push_back(i);
        }

        std::vector<std::vector<std::vector<int>>> partial = {{}};
        for (int pi : part_indices) {
            const auto& part = partition[pi];
            const int n = static_cast<int>(part.size());
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

    return combinations;
}

bool is_subset(const std::vector<int>& small, const std::vector<int>& big) {
    return std::includes(big.begin(), big.end(), small.begin(), small.end());
}

// Barycenter positions per combination: the centroid of the barycenters of
// the parts of the combination's sub-partition. For 2-color partitions this
// equals the flat average of the original points; for 3+ colors it differs.
std::vector<Point3D> combination_barycenters(
    const std::vector<Combination>& combinations,
    const Points& points
) {
    std::vector<Point3D> barycenters;
    barycenters.reserve(combinations.size());
    for (const auto& comb : combinations) {
        Point3D sum = Point3D::Zero();
        for (const auto& part : comb.parts) {
            sum += compute_barycenter(points, part);
        }
        barycenters.push_back(sum / static_cast<double>(comb.parts.size()));
    }
    return barycenters;
}

// Maximal chains of combinations under flat-set inclusion, one element per
// level. Only defined when there are >= 3 distinct levels; each chain becomes
// a higher-dimensional simplex of the subdivision.
std::vector<std::vector<size_t>> enumerate_maximal_chains(
    const std::vector<Combination>& combinations
) {
    std::map<size_t, std::vector<size_t>> levels;
    for (size_t i = 0; i < combinations.size(); ++i) {
        levels[combinations[i].flat.size()].push_back(i);
    }
    if (levels.size() < 3) return {};

    std::vector<std::vector<size_t>> chains;
    auto level_it = levels.begin();
    for (size_t idx : level_it->second) {
        chains.push_back({idx});
    }

    for (++level_it; level_it != levels.end(); ++level_it) {
        std::vector<std::vector<size_t>> new_chains;
        for (const auto& chain : chains) {
            size_t last = chain.back();
            for (size_t idx : level_it->second) {
                if (is_subset(combinations[last].flat, combinations[idx].flat)) {
                    auto extended = chain;
                    extended.push_back(idx);
                    new_chains.push_back(std::move(extended));
                }
            }
        }
        chains = std::move(new_chains);
    }

    return chains;
}

} // namespace

BarycentricSubdivision::BarycentricSubdivision(
    const Points& points,
    const ColorLabels& color_labels,
    bool lower_star
) : points_(points), color_labels_(color_labels), lower_star_(lower_star) {}

double BarycentricSubdivision::star_value(double a, double b) const {
    return lower_star_ ? std::max(a, b) : std::min(a, b);
}

double BarycentricSubdivision::star_value(const std::vector<double>& vals) const {
    if (vals.empty()) return 0.0;
    double result = vals[0];
    for (size_t i = 1; i < vals.size(); ++i) {
        result = star_value(result, vals[i]);
    }
    return result;
}

Point3D BarycentricSubdivision::get_barycenter(const std::vector<int>& vertices) const {
    return compute_barycenter(points_, vertices);
}

double BarycentricSubdivision::compute_filtration_value(const Partition& partitioning) const {
    const size_t k = partitioning.size();
    if (k < 2) return 0.0;

    std::vector<Point3D> bcs;
    bcs.reserve(k);
    for (const auto& part : partitioning) {
        bcs.push_back(get_barycenter(part));
    }

    double sum = 0.0;
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            sum += euclidean_distance(bcs[i], bcs[j]);
        }
    }
    return sum / (k * (k - 1) / 2.0);
}

BarycentricSubdivision::SimplexInfo BarycentricSubdivision::get_or_create_simplex(
    const std::vector<std::vector<int>>& partitioning
) {
    std::vector<int> key;
    for (const auto& part : partitioning) {
        key.insert(key.end(), part.begin(), part.end());
    }
    std::sort(key.begin(), key.end());

    auto it = simplex_map_.lower_bound(key);
    if (it != simplex_map_.end() && it->first == key) {
        return SimplexInfo{it->second.first, it->second.second, false};
    }

    int32_t id = next_simplex_id_++;
    double value = compute_filtration_value(partitioning);
    simplex_map_.emplace_hint(it, std::move(key), std::make_pair(id, value));
    return SimplexInfo{id, value, true};
}

void BarycentricSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    auto partition = get_chromatic_partitioning(simplex_vertices, color_labels_);
    if (partition.size() < 2) return;

    auto combinations = enumerate_multicolored_combinations(partition);
    const size_t n = combinations.size();

    std::vector<std::pair<int32_t, double>> vertices(n);
    std::vector<char> created(n);
    for (size_t i = 0; i < n; ++i) {
        auto info = get_or_create_simplex(combinations[i].parts);
        vertices[i] = {info.id, info.value};
        created[i] = info.newly_created;
    }

    auto new_barycenters = combination_barycenters(combinations, points_);
    for (size_t i = 0; i < n; ++i) {
        if (created[i]) {
            barycenters_.push_back(new_barycenters[i]);
        }
    }

    // Edges: one per strict flat-set inclusion between combinations.
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && is_subset(combinations[i].flat, combinations[j].flat)) {
                Simplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                filtration_set_.insert({edge, star_value(vertices[i].second, vertices[j].second)});
            }
        }
    }

    for (const auto& chain : enumerate_maximal_chains(combinations)) {
        Simplex simplex;
        std::vector<double> vals;
        for (size_t idx : chain) {
            simplex.push_back(vertices[idx].first);
            vals.push_back(vertices[idx].second);
        }
        std::sort(simplex.begin(), simplex.end());
        filtration_set_.insert({simplex, star_value(vals)});
    }

    for (const auto& [id, val] : vertices) {
        filtration_set_.insert({{id}, val});
    }
}

std::vector<std::vector<int>> BarycentricSubdivision::get_vertex_atom_indices() const {
    std::vector<std::vector<int>> result(next_simplex_id_);
    for (const auto& [key, id_val] : simplex_map_) {
        result[id_val.first] = key;
    }
    return result;
}

void BarycentricSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    process_simplex({tet.begin(), tet.end()});
}

Filtration BarycentricSubdivision::get_filtration() const {
    Filtration result(filtration_set_.begin(), filtration_set_.end());

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

std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    bool weighted,
    bool alpha,
    bool lower_star
) {
    detail::validate_inputs(points, color_labels, radii, weighted);

    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    auto mc_simplices = detail::run_subdivision(subdivision, points, color_labels, radii, weighted, alpha);

    return {subdivision.get_barycenters(), subdivision.get_filtration(), std::move(mc_simplices),
            subdivision.get_vertex_atom_indices()};
}

} // namespace delaunay_interfaces
