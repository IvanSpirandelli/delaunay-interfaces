#pragma once

#include "types.hpp"
#include <map>
#include <set>

namespace delaunay_interfaces {

// Barycentric subdivision helper class
class BarycentricSubdivision {
public:
    BarycentricSubdivision(const Points& points, const ColorLabels& color_labels, bool lower_star = false);

    // Process a single tetrahedron (delegates to process_simplex)
    void process_tetrahedron(const Tetrahedron& tet);

    // Process any simplex (generic: works for edges, triangles, tetrahedra)
    void process_simplex(const std::vector<int>& simplex_vertices);

    // Get results
    const Points& get_barycenters() const { return barycenters_; }
    Filtration get_filtration() const;

    // Get atom indices for each barycenter vertex, ordered by vertex ID.
    // vertex_atom_indices[i] = sorted list of input atom indices that vertex i is the barycenter of.
    std::vector<std::vector<int>> get_vertex_atom_indices() const;

private:
    // Barycenter computation
    Point3D get_barycenter(const std::vector<int>& vertices) const;
    Point3D get_barycenter_from_points(const std::vector<Point3D>& points) const;

    // Filtration value computation
    double compute_filtration_value(const Partition& partitioning) const;

    // Get or create barycenter simplex
    struct SimplexInfo {
        int32_t id;
        double value;
        bool newly_created;
    };

    SimplexInfo get_or_create_simplex(const std::vector<std::vector<int>>& partitioning);

    // Star filtration helpers
    double star_value(double a, double b) const;
    double star_value(std::initializer_list<double> vals) const;
    double star_value(const std::vector<double>& vals) const;

    // Data members
    const Points& points_;
    const ColorLabels& color_labels_;
    bool lower_star_;
    Points barycenters_;

    // Map from sorted vertex sets to (simplex_id, filtration_value)
    std::map<std::vector<int>, std::pair<int32_t, double>> simplex_map_;
    int32_t next_simplex_id_ = 0;

    // Filtration simplices
    std::set<SimplexWithFiltration> filtration_set_;
};

} // namespace delaunay_interfaces
