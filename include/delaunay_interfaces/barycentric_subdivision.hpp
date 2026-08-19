#pragma once

#include "types.hpp"
#include <array>
#include <cstdint>
#include <unordered_map>

namespace delaunay_interfaces {

class BarycentricSubdivision {
public:
    // points and color_labels are held by reference and must outlive this object.
    BarycentricSubdivision(const Points& points, const ColorLabels& color_labels, bool lower_star = false);

    void process_tetrahedron(const Tetrahedron& tet);

    // Generic over simplex dimension: works for edges, triangles, and tetrahedra.
    void process_simplex(const std::vector<int>& simplex_vertices);

    // Capacity hint before processing: expected vertex count (an upper bound
    // is fine) and edge / chain filtration entries to be pushed.
    void reserve(size_t vertices, size_t edges, size_t chains) {
        barycenters_.reserve(vertices);
        vertex_map_.reserve(vertices);
        flat_.dim1.reserve(edges);
        flat_.dim2.reserve(chains);
    }

    [[nodiscard]] const Points& get_barycenters() const { return barycenters_; }

    // Sorts and deduplicates the accumulated entries (once), then returns
    // a copy of the lazily materialized filtration; safe to call repeatedly.
    [[nodiscard]] Filtration get_filtration();

    // Move-out variants for result assembly; the moved-from field is spent.
    [[nodiscard]] Points take_barycenters() { return std::move(barycenters_); }
    [[nodiscard]] Filtration take_filtration();

    // Move-out of the sorted flat form without materializing the
    // vector-of-vectors; the flat storage (and any cache) is spent.
    [[nodiscard]] FlatFiltration take_flat_filtration();

    // Atom indices for each barycenter vertex, ordered by vertex ID.
    // vertex_atom_indices[i] = sorted list of input atom indices that vertex i is the barycenter of.
    [[nodiscard]] std::vector<std::vector<int>> get_vertex_atom_indices() const;

private:
    struct VertexInfo {
        int32_t id;
        double value;
    };

    void finalize_filtration();

    // Shared allocation-free core of process_simplex / process_tetrahedron.
    void process_simplex_impl(const int* verts, int n_verts);

    // atoms must be the sorted union of the face's parts (at most n_atoms
    // entries); it doubles as the vertex-map key. parts[p] holds the first
    // part_sizes[p] atom ids of part p, in part order. Creation also appends
    // the vertex's barycenter.
    VertexInfo get_or_create_vertex(
        const int* atoms, int n_atoms,
        const int parts[][4], const int8_t* part_sizes, int n_parts);

    // Upper star (default) takes the minimum over values, lower star the maximum.
    double star_value(double a, double b) const;
    double star_value(const std::vector<double>& vals) const;

    const Points& points_;
    const ColorLabels& color_labels_;
    bool lower_star_;
    Points barycenters_;

    // Sorted atom set packed into a fixed-size key, padded with -1. Simplices
    // have at most 4 vertices (see process_simplex), so atom sets fit.
    using AtomKey = std::array<int, 4>;
    struct AtomKeyHash {
        size_t operator()(const AtomKey& k) const {
            uint64_t h = 0;
            for (int a : k) {
                // splitmix64 round per element
                uint64_t x = h ^ (static_cast<uint64_t>(static_cast<uint32_t>(a)) + 0x9e3779b97f4a7c15ULL);
                x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
                x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
                h = x ^ (x >> 31);
            }
            return static_cast<size_t>(h);
        }
    };

    // Map from packed sorted atom sets to (vertex_id, filtration_value)
    std::unordered_map<AtomKey, std::pair<int32_t, double>, AtomKeyHash> vertex_map_;
    int32_t next_vertex_id_ = 0;

    // Per-dimension POD accumulation buckets (every emitted simplex has 1, 2,
    // or 3 vertices), stored directly in the FlatFiltration that becomes the
    // surface's primary filtration form. flat_.dim1/dim2 may hold duplicates
    // from shared faces of adjacent tetrahedra; duplicate entries are exact
    // copies (values derive from the vertex ids alone), so
    // finalize_filtration's per-bucket sort + unique removes them. flat_.dim0
    // (vertex singletons) is filled from vertex_map_ during finalization and
    // is duplicate-free by construction.
    FlatFiltration flat_;
    bool finalized_ = false;
};

} // namespace delaunay_interfaces
