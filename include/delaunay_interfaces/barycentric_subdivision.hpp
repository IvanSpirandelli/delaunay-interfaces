#pragma once

#include "types.hpp"
#include <map>

namespace delaunay_interfaces {

class BarycentricSubdivision {
public:
    // points and color_labels are held by reference and must outlive this object.
    BarycentricSubdivision(const Points& points, const ColorLabels& color_labels, bool lower_star = false);

    void process_tetrahedron(const Tetrahedron& tet);

    // Generic over simplex dimension: works for edges, triangles, and tetrahedra.
    void process_simplex(const std::vector<int>& simplex_vertices);

    [[nodiscard]] const Points& get_barycenters() const { return barycenters_; }

    // Sorts and deduplicates the accumulated entries in place, then returns
    // a copy; safe to call repeatedly.
    [[nodiscard]] Filtration get_filtration();

    // Move-out variants for result assembly; the moved-from field is spent.
    [[nodiscard]] Points take_barycenters() { return std::move(barycenters_); }
    [[nodiscard]] Filtration take_filtration();

    // Atom indices for each barycenter vertex, ordered by vertex ID.
    // vertex_atom_indices[i] = sorted list of input atom indices that vertex i is the barycenter of.
    [[nodiscard]] std::vector<std::vector<int>> get_vertex_atom_indices() const;

private:
    struct VertexInfo {
        int32_t id;
        double value;
    };

    Point3D get_barycenter(const std::vector<int>& vertices) const;
    void finalize_filtration();

    // atoms must be the sorted union of the partition's parts; it doubles as
    // the vertex-map key. Creation also appends the vertex's barycenter.
    VertexInfo get_or_create_vertex(const Partition& partition, const std::vector<int>& atoms);

    // Upper star (default) takes the minimum over values, lower star the maximum.
    double star_value(double a, double b) const;
    double star_value(const std::vector<double>& vals) const;

    const Points& points_;
    const ColorLabels& color_labels_;
    bool lower_star_;
    Points barycenters_;

    // Map from sorted atom sets to (vertex_id, filtration_value)
    std::map<std::vector<int>, std::pair<int32_t, double>> vertex_map_;
    int32_t next_vertex_id_ = 0;

    // May hold duplicates from shared faces of adjacent tetrahedra; duplicate
    // entries are exact copies (values derive from the vertex ids alone), so
    // get_filtration's sort + unique removes them.
    Filtration filtration_;
};

} // namespace delaunay_interfaces
