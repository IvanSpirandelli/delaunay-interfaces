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
                int32_t a = vertices[i].first;
                int32_t b = vertices[j].first;
                if (a > b) std::swap(a, b);
                dim2_.push_back({star_value(vertices[i].second, vertices[j].second), {a, b}});
            }
        }
    }

    // Chains always have exactly 3 elements: one per distinct atom-set size,
    // and the sizes lie in {2, 3, 4} (>= 2 atoms per face, <= 4 in total).
    for (const auto& chain : enumerate_maximal_chains(faces)) {
        int32_t ids[3] = {
            vertices[chain[0]].first, vertices[chain[1]].first, vertices[chain[2]].first};
        double value = star_value(
            star_value(vertices[chain[0]].second, vertices[chain[1]].second),
            vertices[chain[2]].second);
        if (ids[0] > ids[1]) std::swap(ids[0], ids[1]);
        if (ids[1] > ids[2]) std::swap(ids[1], ids[2]);
        if (ids[0] > ids[1]) std::swap(ids[0], ids[1]);
        dim3_.push_back({value, {ids[0], ids[1], ids[2]}});
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

// Sort each bucket by (value, simplex) — a total order, making the output
// deterministic and equal entries adjacent so unique can drop the duplicates
// from shared faces — then materialize the public Filtration once as
// dim1|dim2|dim3, which reproduces the global (dimension, value, simplex)
// order of the former single-vector sort exactly.
void BarycentricSubdivision::finalize_filtration() {
    if (finalized_) return;
    finalized_ = true;

    // Vertex singletons are emitted here, once, instead of per processed
    // simplex: every map entry was created for some simplex that would have
    // pushed its singleton, so the entry set (and hence the sorted, deduped
    // output) is unchanged. Duplicate-free by construction, so no unique.
    dim1_.reserve(vertex_map_.size());
    for (const auto& [key, id_val] : vertex_map_) {
        (void)key;
        dim1_.push_back({id_val.second, id_val.first});
    }
    std::sort(dim1_.begin(), dim1_.end(),
        [](const Dim1Entry& a, const Dim1Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            return a.v < b.v;
        });

    std::sort(dim2_.begin(), dim2_.end(),
        [](const Dim2Entry& a, const Dim2Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            return a.v[1] < b.v[1];
        });
    dim2_.erase(std::unique(dim2_.begin(), dim2_.end(),
        [](const Dim2Entry& a, const Dim2Entry& b) {
            return a.value == b.value && a.v[0] == b.v[0] && a.v[1] == b.v[1];
        }), dim2_.end());

    std::sort(dim3_.begin(), dim3_.end(),
        [](const Dim3Entry& a, const Dim3Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            if (a.v[1] != b.v[1]) return a.v[1] < b.v[1];
            return a.v[2] < b.v[2];
        });
    dim3_.erase(std::unique(dim3_.begin(), dim3_.end(),
        [](const Dim3Entry& a, const Dim3Entry& b) {
            return a.value == b.value && a.v[0] == b.v[0] && a.v[1] == b.v[1]
                && a.v[2] == b.v[2];
        }), dim3_.end());

    filtration_.reserve(dim1_.size() + dim2_.size() + dim3_.size());
    for (const auto& e : dim1_) {
        filtration_.emplace_back(SurfaceSimplex{e.v}, e.value);
    }
    for (const auto& e : dim2_) {
        filtration_.emplace_back(SurfaceSimplex{e.v[0], e.v[1]}, e.value);
    }
    for (const auto& e : dim3_) {
        filtration_.emplace_back(SurfaceSimplex{e.v[0], e.v[1], e.v[2]}, e.value);
    }

    // The buckets are spent; release their memory.
    dim1_ = {};
    dim2_ = {};
    dim3_ = {};
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
