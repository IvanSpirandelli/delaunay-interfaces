#include "delaunay_interfaces/barycentric_subdivision.hpp"
#include "delaunay_interfaces/chromatic_partitioning.hpp"
#include "subdivision_driver.hpp"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
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

// Only 7 partition shapes (ordered part-size sequences as
// compute_chromatic_partition orders them) are reachable from multicolored
// simplices with <= 4 vertices: tets (3,1), (2,2), (2,1,1), (1,1,1,1);
// triangles (2,1), (1,1,1); edges (1,1). The face enumeration, the inclusion
// pairs, and the maximal chains depend only on the shape, never on the actual
// atom ids or colors, so they are precomputed once per shape by running the
// enumerators above on a canonical representative whose atoms are the
// flattened part positions 0..n-1. Replaying the stored enumeration order at
// runtime keeps vertex creation order, barycenter summation order, star-value
// fold order, and emission order identical to the direct enumeration, so the
// output stays bitwise identical.
constexpr int kNumShapes = 7;
constexpr int kMaxFaces = 11; // shape (1,1,1,1)

struct FaceEntry {
    int8_t n_parts;
    int8_t part_sizes[4];
    int8_t part_positions[4][4]; // canonical positions per part, in part order
    int8_t n_atoms;
    int8_t atom_positions[4];    // sorted union of the part positions
};

struct ShapeTable {
    std::vector<FaceEntry> faces;                   // in enumeration order
    std::vector<std::array<int8_t, 2>> inclusion_pairs; // (i, j) face indices
    std::vector<std::array<int8_t, 3>> chains;      // face-index triples
};

// Index must match the shape list in build_shape_tables.
int shape_index(int n_parts, const int* sizes) {
    switch (n_parts) {
        case 2:
            if (sizes[0] == 1) return 0;                     // (1,1)
            if (sizes[0] == 2) return sizes[1] == 1 ? 1 : 4; // (2,1) / (2,2)
            return 3;                                        // (3,1)
        case 3:
            return sizes[0] == 2 ? 5 : 2;                    // (2,1,1) / (1,1,1)
        default:
            return 6;                                        // (1,1,1,1)
    }
}

std::array<ShapeTable, kNumShapes> build_shape_tables() {
    const std::vector<std::vector<int>> shapes = {
        {1, 1}, {2, 1}, {1, 1, 1}, {3, 1}, {2, 2}, {2, 1, 1}, {1, 1, 1, 1}};

    std::array<ShapeTable, kNumShapes> tables;
    for (int s = 0; s < kNumShapes; ++s) {
        Partition partition;
        int pos = 0;
        for (int size : shapes[s]) {
            std::vector<int> part;
            for (int i = 0; i < size; ++i) part.push_back(pos++);
            partition.push_back(std::move(part));
        }

        auto faces = enumerate_multicolored_faces(partition);
        ShapeTable& table = tables[s];
        if (faces.size() > static_cast<size_t>(kMaxFaces)) {
            throw std::logic_error("shape table exceeds kMaxFaces");
        }
        for (const auto& face : faces) {
            FaceEntry e{};
            e.n_parts = static_cast<int8_t>(face.parts.size());
            for (size_t p = 0; p < face.parts.size(); ++p) {
                e.part_sizes[p] = static_cast<int8_t>(face.parts[p].size());
                for (size_t i = 0; i < face.parts[p].size(); ++i) {
                    e.part_positions[p][i] = static_cast<int8_t>(face.parts[p][i]);
                }
            }
            e.n_atoms = static_cast<int8_t>(face.atoms.size());
            for (size_t i = 0; i < face.atoms.size(); ++i) {
                e.atom_positions[i] = static_cast<int8_t>(face.atoms[i]);
            }
            table.faces.push_back(e);
        }

        // Strict-inclusion pairs, in the exact (i outer, j inner) order the
        // former per-simplex double loop found them.
        const size_t n = faces.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i != j && is_subset(faces[i].atoms, faces[j].atoms)) {
                    table.inclusion_pairs.push_back(
                        {static_cast<int8_t>(i), static_cast<int8_t>(j)});
                }
            }
        }

        for (const auto& chain : enumerate_maximal_chains(faces)) {
            if (chain.size() != 3) {
                throw std::logic_error("maximal chain without exactly 3 levels");
            }
            table.chains.push_back({static_cast<int8_t>(chain[0]),
                static_cast<int8_t>(chain[1]), static_cast<int8_t>(chain[2])});
        }
    }
    return tables;
}

const std::array<ShapeTable, kNumShapes>& shape_tables() {
    static const std::array<ShapeTable, kNumShapes> tables = build_shape_tables();
    return tables;
}

// Sorts a filtration bucket by (value, ids lexicographic), matching
// std::sort with the same comparator exactly. Filtration values are mean
// pairwise distances: always >= 0 with a clear sign bit and never NaN, so
// the raw 64-bit IEEE-754 pattern ordered as an unsigned integer is exactly
// the numeric order, and a stable LSD radix over its 8 bytes sorts by value.
// The radix leaves equal-value entries in emission order, not id order, so
// each equal-value run is re-sorted by ids afterwards (equal values here
// means equal bit patterns — duplicates derive from identical arithmetic —
// and the runs are typically tiny). A pass whose byte is shared by every key
// is the identity permutation and is skipped; small buckets fall back to
// std::sort, which is faster there and keeps tiny cases fast.
template <typename Entry, typename FullLess, typename IdLess>
void sort_bucket(std::vector<Entry>& entries, FullLess full_less, IdLess id_less) {
    constexpr size_t kRadixMinSize = 10000;
    const size_t n = entries.size();
    if (n < kRadixMinSize) {
        std::sort(entries.begin(), entries.end(), full_less);
        return;
    }

#ifndef NDEBUG
    for (const Entry& e : entries) {
        assert(e.value >= 0.0 && !std::signbit(e.value));
    }
#endif

    auto key_bits = [](const Entry& e) {
        uint64_t bits;
        std::memcpy(&bits, &e.value, sizeof bits);
        return bits;
    };

    std::vector<Entry> scratch(n);
    Entry* src = entries.data();
    Entry* dst = scratch.data();
    for (int shift = 0; shift < 64; shift += 8) {
        size_t count[256] = {};
        for (size_t i = 0; i < n; ++i) {
            ++count[(key_bits(src[i]) >> shift) & 0xff];
        }
        if (std::find(std::begin(count), std::end(count), n) != std::end(count)) {
            continue;
        }
        size_t offset = 0;
        for (size_t& c : count) {
            const size_t bucket = c;
            c = offset;
            offset += bucket;
        }
        for (size_t i = 0; i < n; ++i) {
            dst[count[(key_bits(src[i]) >> shift) & 0xff]++] = src[i];
        }
        std::swap(src, dst);
    }
    if (src != entries.data()) {
        std::memcpy(entries.data(), src, n * sizeof(Entry));
    }

    size_t run_begin = 0;
    for (size_t i = 1; i <= n; ++i) {
        if (i == n || entries[i].value != entries[run_begin].value) {
            if (i - run_begin > 1) {
                std::sort(entries.begin() + run_begin, entries.begin() + i, id_less);
            }
            run_begin = i;
        }
    }
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

// Barycenter position per combination: the centroid of the barycenters of
// the parts of the combination's sub-partition. This equals the flat average
// of the original points only when all chosen parts have equal size (e.g.
// 1-1, 2-2); for unequal part sizes (2-1, 3-1, ...) it differs.
BarycentricSubdivision::VertexInfo BarycentricSubdivision::get_or_create_vertex(
    const int* atoms, int n_atoms,
    const int parts[][4], const int8_t* part_sizes, int n_parts
) {
    AtomKey key;
    key.fill(-1);
    std::copy(atoms, atoms + n_atoms, key.begin());
    auto it = vertex_map_.find(key);
    if (it != vertex_map_.end()) {
        return VertexInfo{it->second.first, it->second.second};
    }

    // The part barycenters yield both the filtration value (mean pairwise
    // distance; faces always choose >= 2 parts) and the vertex position.
    // Summation orders match the former compute_barycenter / pairwise loops
    // exactly, keeping the results bitwise identical.
    const size_t k = static_cast<size_t>(n_parts);
    Point3D bcs[4];
    for (int p = 0; p < n_parts; ++p) {
        Point3D center = Point3D::Zero();
        for (int s = 0; s < part_sizes[p]; ++s) {
            center += points_[parts[p][s]];
        }
        bcs[p] = center / static_cast<double>(part_sizes[p]);
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
    process_simplex_impl(simplex_vertices.data(), static_cast<int>(simplex_vertices.size()));
}

void BarycentricSubdivision::process_simplex_impl(const int* verts, int n_verts) {
    if (n_verts > 4) {
        throw std::invalid_argument(
            "process_simplex supports at most 4 vertices (edge, triangle, or tetrahedron)");
    }

    // Reproduce compute_chromatic_partition on fixed arrays: parts grouped in
    // ascending color order (the std::map iteration order) with atoms appended
    // in vertex order, then a stable sort of the parts by descending size.
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
    if (n_parts < 2) return;

    // Stable insertion sort by descending size (ties keep color order).
    int order[4] = {0, 1, 2, 3};
    for (int i = 1; i < n_parts; ++i) {
        const int key = order[i];
        int j = i;
        while (j > 0 && part_sizes[order[j - 1]] < part_sizes[key]) {
            order[j] = order[j - 1];
            --j;
        }
        order[j] = key;
    }

    // Flattened position -> actual atom id, in part order; the canonical
    // shape representatives use these positions as their atoms.
    int sizes[4];
    int pos_to_id[4];
    int flat = 0;
    for (int i = 0; i < n_parts; ++i) {
        const int p = order[i];
        sizes[i] = part_sizes[p];
        for (int s = 0; s < part_sizes[p]; ++s) {
            pos_to_id[flat++] = part_atoms[p][s];
        }
    }

    const ShapeTable& table = shape_tables()[shape_index(n_parts, sizes)];
    const int n = static_cast<int>(table.faces.size());

    std::pair<int32_t, double> vertices[kMaxFaces];
    for (int i = 0; i < n; ++i) {
        const FaceEntry& face = table.faces[i];
        // pos_to_id is not monotonic, so the atom ids must be re-sorted to
        // match the sorted unions the direct enumeration produced.
        int atoms[4];
        for (int a = 0; a < face.n_atoms; ++a) {
            atoms[a] = pos_to_id[face.atom_positions[a]];
        }
        for (int a = 1; a < face.n_atoms; ++a) {
            const int key = atoms[a];
            int b = a;
            while (b > 0 && atoms[b - 1] > key) {
                atoms[b] = atoms[b - 1];
                --b;
            }
            atoms[b] = key;
        }
        int face_parts[4][4];
        for (int p = 0; p < face.n_parts; ++p) {
            for (int s = 0; s < face.part_sizes[p]; ++s) {
                face_parts[p][s] = pos_to_id[face.part_positions[p][s]];
            }
        }
        auto info = get_or_create_vertex(
            atoms, face.n_atoms, face_parts, face.part_sizes, face.n_parts);
        vertices[i] = {info.id, info.value};
    }

    // Edges: one per strict inclusion between the faces' atom sets.
    for (const auto& pair : table.inclusion_pairs) {
        const auto& vi = vertices[pair[0]];
        const auto& vj = vertices[pair[1]];
        int32_t a = vi.first;
        int32_t b = vj.first;
        if (a > b) std::swap(a, b);
        flat_.dim2.push_back({star_value(vi.second, vj.second), {a, b}});
    }

    // Chains always have exactly 3 elements: one per distinct atom-set size,
    // and the sizes lie in {2, 3, 4} (>= 2 atoms per face, <= 4 in total).
    for (const auto& chain : table.chains) {
        int32_t ids[3] = {
            vertices[chain[0]].first, vertices[chain[1]].first, vertices[chain[2]].first};
        double value = star_value(
            star_value(vertices[chain[0]].second, vertices[chain[1]].second),
            vertices[chain[2]].second);
        if (ids[0] > ids[1]) std::swap(ids[0], ids[1]);
        if (ids[1] > ids[2]) std::swap(ids[1], ids[2]);
        if (ids[0] > ids[1]) std::swap(ids[0], ids[1]);
        flat_.dim3.push_back({value, {ids[0], ids[1], ids[2]}});
    }
}

void BarycentricSubdivision::process_tetrahedron(const Tetrahedron& tet) {
    process_simplex_impl(tet.data(), 4);
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
// from shared faces. The sorted buckets ARE the result (flat_ is the primary
// form); the public vector-of-vectors Filtration is materialized lazily by
// FlatFiltration as dim1|dim2|dim3, which reproduces the global
// (dimension, value, simplex) order of the former single-vector sort exactly.
void BarycentricSubdivision::finalize_filtration() {
    if (finalized_) return;
    finalized_ = true;

    using Dim1Entry = FlatFiltration::Dim1Entry;
    using Dim2Entry = FlatFiltration::Dim2Entry;
    using Dim3Entry = FlatFiltration::Dim3Entry;

    // Vertex singletons are emitted here, once, instead of per processed
    // simplex: every map entry was created for some simplex that would have
    // pushed its singleton, so the entry set (and hence the sorted, deduped
    // output) is unchanged. Duplicate-free by construction, so no unique.
    flat_.dim1.reserve(vertex_map_.size());
    for (const auto& [key, id_val] : vertex_map_) {
        (void)key;
        flat_.dim1.push_back({id_val.second, id_val.first});
    }
    sort_bucket(flat_.dim1,
        [](const Dim1Entry& a, const Dim1Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            return a.v < b.v;
        },
        [](const Dim1Entry& a, const Dim1Entry& b) { return a.v < b.v; });

    sort_bucket(flat_.dim2,
        [](const Dim2Entry& a, const Dim2Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            return a.v[1] < b.v[1];
        },
        [](const Dim2Entry& a, const Dim2Entry& b) {
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            return a.v[1] < b.v[1];
        });
    flat_.dim2.erase(std::unique(flat_.dim2.begin(), flat_.dim2.end(),
        [](const Dim2Entry& a, const Dim2Entry& b) {
            return a.value == b.value && a.v[0] == b.v[0] && a.v[1] == b.v[1];
        }), flat_.dim2.end());

    sort_bucket(flat_.dim3,
        [](const Dim3Entry& a, const Dim3Entry& b) {
            if (a.value != b.value) return a.value < b.value;
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            if (a.v[1] != b.v[1]) return a.v[1] < b.v[1];
            return a.v[2] < b.v[2];
        },
        [](const Dim3Entry& a, const Dim3Entry& b) {
            if (a.v[0] != b.v[0]) return a.v[0] < b.v[0];
            if (a.v[1] != b.v[1]) return a.v[1] < b.v[1];
            return a.v[2] < b.v[2];
        });
    flat_.dim3.erase(std::unique(flat_.dim3.begin(), flat_.dim3.end(),
        [](const Dim3Entry& a, const Dim3Entry& b) {
            return a.value == b.value && a.v[0] == b.v[0] && a.v[1] == b.v[1]
                && a.v[2] == b.v[2];
        }), flat_.dim3.end());
}

Filtration BarycentricSubdivision::get_filtration() {
    finalize_filtration();
    return flat_.materialized();
}

Filtration BarycentricSubdivision::take_filtration() {
    finalize_filtration();
    return flat_.take_materialized();
}

FlatFiltration BarycentricSubdivision::take_flat_filtration() {
    finalize_filtration();
    return std::move(flat_);
}

namespace detail {

std::tuple<Points, FlatFiltration, MulticoloredSimplices, VertexAtomIndices> compute_subdivision_flat(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha,
    bool lower_star
) {
    validate_inputs(points, color_labels, radii);

    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    auto mc_simplices = run_subdivision(
        subdivision, points, color_labels, radii, alpha.value_or(!radii.empty()));

    return {subdivision.take_barycenters(), subdivision.take_flat_filtration(), std::move(mc_simplices),
            subdivision.get_vertex_atom_indices()};
}

std::tuple<Points, FlatFiltration, MulticoloredSimplices, VertexAtomIndices> compute_subdivision_flat(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool lower_star
) {
    validate_inputs(points, color_labels, Radii{});

    BarycentricSubdivision subdivision(points, color_labels, lower_star);
    auto mc_simplices = run_subdivision(subdivision, points, color_labels, radius);

    return {subdivision.take_barycenters(), subdivision.take_flat_filtration(), std::move(mc_simplices),
            subdivision.get_vertex_atom_indices()};
}

} // namespace detail

std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> alpha,
    bool lower_star
) {
    auto [vertices, flat, mc_simplices, vertex_atom_indices] = detail::compute_subdivision_flat(
        points, color_labels, radii, alpha, lower_star);
    return {std::move(vertices), flat.take_materialized(), std::move(mc_simplices),
            std::move(vertex_atom_indices)};
}

std::tuple<Points, Filtration, MulticoloredSimplices, VertexAtomIndices> compute_barycentric_subdivision_and_filtration(
    const Points& points,
    const ColorLabels& color_labels,
    double radius,
    bool lower_star
) {
    auto [vertices, flat, mc_simplices, vertex_atom_indices] = detail::compute_subdivision_flat(
        points, color_labels, radius, lower_star);
    return {std::move(vertices), flat.take_materialized(), std::move(mc_simplices),
            std::move(vertex_atom_indices)};
}

} // namespace delaunay_interfaces
