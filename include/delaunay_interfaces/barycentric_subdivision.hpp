#pragma once

#include "types.hpp"
#include <map>
#include <set>

namespace delaunay_interfaces {

// Barycentric subdivision helper class
class BarycentricSubdivision {
public:
    BarycentricSubdivision(const Points& points, const ColorLabels& color_labels);

    // Process a single tetrahedron
    void process_tetrahedron(const Tetrahedron& tet);

    // Get results
    const Points& get_barycenters() const { return barycenters_; }
    Filtration get_filtration() const;

private:
    // Chromatic partitioning
    Partition get_chromatic_partitioning(const Tetrahedron& tet) const;

    // Barycenter computation
    Point3D get_barycenter(const std::vector<int>& vertices) const;
    Point3D get_barycenter_from_points(const std::vector<Point3D>& points) const;

    // Filtration value computation
    double compute_filtration_value(const Partition& partitioning) const;

    // Scaffold extension for different partition types
    void extend_scaffold_2_2(const std::vector<int>& part1, const std::vector<int>& part2);
    void extend_scaffold_3_1(const std::vector<int>& part1, const std::vector<int>& part2);
    void extend_scaffold_2_1_1(
        const std::vector<int>& part1,
        const std::vector<int>& part2,
        const std::vector<int>& part3
    );
    void extend_scaffold_1_1_1_1(
        const std::vector<int>& part1,
        const std::vector<int>& part2,
        const std::vector<int>& part3,
        const std::vector<int>& part4
    );

    // Get or create barycenter simplex
    struct SimplexInfo {
        int32_t id;
        double value;
        bool newly_created;
    };

    SimplexInfo get_or_create_simplex(const std::vector<std::vector<int>>& partitioning);

    // Data members
    const Points& points_;
    const ColorLabels& color_labels_;
    Points barycenters_;

    // Map from sorted vertex sets to (simplex_id, filtration_value)
    std::map<std::vector<int>, std::pair<int32_t, double>> simplex_map_;
    int32_t next_simplex_id_ = 0;

    // Filtration simplices
    std::set<SimplexWithFiltration> filtration_set_;
};

} // namespace delaunay_interfaces
