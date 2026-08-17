#include "delaunay_interfaces/barycentric_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "subdivision_driver.hpp"
#include <algorithm>
#include <bitset>
#include <stdexcept>

namespace delaunay_interfaces {

namespace {

// One multicolored face of the processed simplex: a choice of >= 2 parts
// and a non-empty subset of each chosen part.
struct MulticoloredFace {
    std::vector<std::vector<int>> parts; // the sub-partition, one entry per chosen part
    std::vector<int> atoms;              // sorted union of the chosen atom indices
};

// Enumerates every subset I of the parts with |I| >= 2, crossed with every
// tuple of non-empty subsets of the parts in I.
std::vector<MulticoloredFace> enumerate_multicolored_faces(const Partition& partition) {
    const int k = static_cast<int>(partition.size());
    std::vector<MulticoloredFace> faces;

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
            MulticoloredFace face;
            face.parts = std::move(tuple);
            for (const auto& s : face.parts) {
                face.atoms.insert(face.atoms.end(), s.begin(), s.end());
            }
            std::sort(face.atoms.begin(), face.atoms.end());
            faces.push_back(std::move(face));
        }
    }

    return faces;
}

bool is_subset(const std::vector<int>& small, const std::vector<int>& big) {
    return std::includes(big.begin(), big.end(), small.begin(), small.end());
}

// Maximal chains of faces under flat-set inclusion, one element per
// level. Only defined when there are >= 3 distinct levels; each chain becomes
// a higher-dimensional simplex of the subdivision.
std::vector<std::vector<size_t>> enumerate_maximal_chains(
    const std::vector<MulticoloredFace>& faces
) {
    std::map<size_t, std::vector<size_t>> levels;
    for (size_t i = 0; i < faces.size(); ++i) {
        levels[faces[i].atoms.size()].push_back(i);
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
                if (is_subset(faces[last].atoms, faces[idx].atoms)) {
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

// Barycenter position per combination: the centroid of the barycenters of
// the parts of the combination's sub-partition. This equals the flat average
// of the original points only when all chosen parts have equal size (e.g.
// 1-1, 2-2); for unequal part sizes (2-1, 3-1, ...) it differs.
BarycentricSubdivision::VertexInfo BarycentricSubdivision::get_or_create_vertex(
    const Partition& partition,
    const std::vector<int>& atoms
) {
    AtomKey key;
    key.fill(-1);
    std::copy(atoms.begin(), atoms.end(), key.begin());
    auto it = vertex_map_.find(key);
    if (it != vertex_map_.end()) {
        return VertexInfo{it->second.first, it->second.second};
    }

    // The part barycenters yield both the filtration value (mean pairwise
    // distance; faces always choose >= 2 parts) and the vertex position.
    const size_t k = partition.size();
    std::vector<Point3D> bcs;
    bcs.reserve(k);
    for (const auto& part : partition) {
        bcs.push_back(get_barycenter(part));
    }

    double sum = 0.0;
    Point3D centroid = Point3D::Zero();
    for (size_t i = 0; i < k; ++i) {
        centroid += bcs[i];
        for (size_t j = i + 1; j < k; ++j) {
            sum += euclidean_distance(bcs[i], bcs[j]);
        }
    }
    double value = sum / (k * (k - 1) / 2.0);
    barycenters_.push_back(centroid / static_cast<double>(k));

    int32_t id = next_vertex_id_++;
    vertex_map_.emplace(key, std::make_pair(id, value));
    return VertexInfo{id, value};
}

void BarycentricSubdivision::process_simplex(const std::vector<int>& simplex_vertices) {
    if (simplex_vertices.size() > 4) {
        throw std::invalid_argument(
            "process_simplex supports at most 4 vertices (edge, triangle, or tetrahedron)");
    }
    auto partition = compute_chromatic_partition(simplex_vertices, color_labels_);
    if (partition.size() < 2) return;

    auto faces = enumerate_multicolored_faces(partition);
    const size_t n = faces.size();

    std::vector<std::pair<int32_t, double>> vertices(n);
    for (size_t i = 0; i < n; ++i) {
        auto info = get_or_create_vertex(faces[i].parts, faces[i].atoms);
        vertices[i] = {info.id, info.value};
    }

    // Edges: one per strict inclusion between the faces' atom sets.
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && is_subset(faces[i].atoms, faces[j].atoms)) {
                SurfaceSimplex edge = {vertices[i].first, vertices[j].first};
                std::sort(edge.begin(), edge.end());
                filtration_.emplace_back(std::move(edge), star_value(vertices[i].second, vertices[j].second));
            }
        }
    }

    for (const auto& chain : enumerate_maximal_chains(faces)) {
        SurfaceSimplex simplex;
        std::vector<double> vals;
        for (size_t idx : chain) {
            simplex.push_back(vertices[idx].first);
            vals.push_back(vertices[idx].second);
        }
        std::sort(simplex.begin(), simplex.end());
        filtration_.emplace_back(std::move(simplex), star_value(vals));
    }
}

void BarycentricSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    process_simplex({tet.begin(), tet.end()});
}

std::vector<std::vector<int>> BarycentricSubdivision::get_vertex_atom_indices() const {
    std::vector<std::vector<int>> result(next_vertex_id_);
    for (const auto& [key, id_val] : vertex_map_) {
        auto& atoms = result[id_val.first];
        for (int a : key) {
            if (a < 0) break;
            atoms.push_back(a);
        }
    }
    return result;
}

// Sort by (dimension, value) with the simplex itself as a total-order
// tiebreak, making the output deterministic and equal entries adjacent
// so unique can drop the duplicates from shared faces.
void BarycentricSubdivision::finalize_filtration() {
    // Vertex singletons are emitted here, once, instead of per processed
    // simplex: every map entry was created for some simplex that would have
    // pushed its singleton, so the entry set (and hence the sorted, deduped
    // output) is unchanged.
    if (!singletons_emitted_) {
        filtration_.reserve(filtration_.size() + vertex_map_.size());
        for (const auto& [key, id_val] : vertex_map_) {
            (void)key;
            filtration_.emplace_back(SurfaceSimplex{id_val.first}, id_val.second);
        }
        singletons_emitted_ = true;
    }
    std::sort(filtration_.begin(), filtration_.end(),
        [](const auto& a, const auto& b) {
            if (std::get<0>(a).size() != std::get<0>(b).size()) {
                return std::get<0>(a).size() < std::get<0>(b).size();
            }
            if (std::get<1>(a) != std::get<1>(b)) {
                return std::get<1>(a) < std::get<1>(b);
            }
            return std::get<0>(a) < std::get<0>(b);
        }
    );
    filtration_.erase(std::unique(filtration_.begin(), filtration_.end()), filtration_.end());
}

Filtration BarycentricSubdivision::get_filtration() {
    finalize_filtration();
    return filtration_;
}

Filtration BarycentricSubdivision::take_filtration() {
    finalize_filtration();
    return std::move(filtration_);
}

std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha,
    bool lower_star
) {
    detail::validate_inputs(points, color_labels, radii);

    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    auto mc_simplices = detail::run_subdivision(
        subdivision, points, color_labels, radii, alpha.value_or(!radii.empty()));

    return {subdivision.take_barycenters(), subdivision.take_filtration(), std::move(mc_simplices),
            subdivision.get_vertex_atom_indices()};
}

std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool lower_star
) {
    detail::validate_inputs(points, color_labels, Radii{});

    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    auto mc_simplices = detail::run_subdivision(subdivision, points, color_labels, radius);

    return {subdivision.take_barycenters(), subdivision.take_filtration(), std::move(mc_simplices),
            subdivision.get_vertex_atom_indices()};
}

} // namespace delaunay_interfaces
