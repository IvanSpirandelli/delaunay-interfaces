#pragma once

#include "types.hpp"
#include <map>
#include <set>

namespace delaunay_interfaces {

class BarycentricSubdivision {
public:
    // points and color_labels are held by reference and must outlive this object.
    BarycentricSubdivision(const Points& points, const ColorLabels& color_labels, bool lower_star = false);

    void process_tetrahedron(const Tetrahedron& tet);

    // Generic over simplex dimension: works for edges, triangles, and tetrahedra.
    void process_simplex(const std::vector<int>& simplex_vertices);

    [[nodiscard]] const Points& get_barycenters() const { return barycenters_; }
    [[nodiscard]] Filtration get_filtration() const;

    // Atom indices for each barycenter vertex, ordered by vertex ID.
    // vertex_atom_indices[i] = sorted list of input atom indices that vertex i is the barycenter of.
    [[nodiscard]] std::vector<std::vector<int>> get_vertex_atom_indices() const;

private:
    struct VertexInfo {
        int32_t id;
        double value;
        bool newly_created;
    };

    Point3D get_barycenter(const std::vector<int>& vertices) const;
    double compute_filtration_value(const Partition& partition) const;
    VertexInfo get_or_create_vertex(const Partition& partition);

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

    std::set<FiltrationEntry> filtration_set_;
};

} // namespace delaunay_interfaces
