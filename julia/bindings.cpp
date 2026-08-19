#include "jlcxx/jlcxx.hpp"
#include "delaunay_interfaces/interface_generation.hpp"
#include <vector>

using namespace delaunay_interfaces;

// Non-mirrored: these types have non-trivial layout, so CxxWrap must box
// them instead of mapping them to isbits Julia structs.
namespace jlcxx {
    template<> struct IsMirroredType<InterfaceGenerator> : std::false_type { };
    template<> struct IsMirroredType<InterfaceSurface> : std::false_type { };
}

// Helper to convert Julia arrays to C++ vectors
template<typename T>
std::vector<T> julia_array_to_vector(jlcxx::ArrayRef<T> arr) {
    return std::vector<T>(arr.begin(), arr.end());
}

// Helper to convert flat point array to C++ Points (points stored as [x1,y1,z1,x2,y2,z2,...])
Points flat_points_to_cpp(jlcxx::ArrayRef<double> flat_points, int n_points) {
    Points points;
    points.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        points.emplace_back(flat_points[i*3], flat_points[i*3+1], flat_points[i*3+2]);
    }
    return points;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    // InterfaceSurface type
    mod.add_type<InterfaceSurface>("InterfaceSurfaceCxx")
        .method("num_vertices", [](const InterfaceSurface& s) {
            return static_cast<int>(s.vertices.size());
        })
        .method("num_simplices", [](const InterfaceSurface& s) {
            return static_cast<int>(s.filtration.size());
        })
        .method("is_alpha", [](const InterfaceSurface& s) {
            return s.alpha;
        })
        .method("is_lower_star", [](const InterfaceSurface& s) {
            return s.lower_star;
        })
        .method("get_vertex", [](const InterfaceSurface& s, int i) {
            if (i < 0 || i >= static_cast<int>(s.vertices.size())) {
                throw std::out_of_range("Vertex index out of range");
            }
            const auto& v = s.vertices[i];
            return std::vector<double>{v[0], v[1], v[2]};
        })
        .method("get_all_vertices", [](const InterfaceSurface& s) {
            std::vector<double> result;
            result.reserve(s.vertices.size() * 3);
            for (const auto& v : s.vertices) {
                result.push_back(v[0]);
                result.push_back(v[1]);
                result.push_back(v[2]);
            }
            return result;
        })
        .method("get_simplex_vertices", [](const InterfaceSurface& s, int i) {
            const auto& f = s.filtration;
            if (i < 0 || i >= static_cast<int>(f.size())) {
                throw std::out_of_range("Simplex index out of range");
            }
            // Entry order is dim0|dim1|dim2, matching the materialized form.
            size_t idx = static_cast<size_t>(i);
            if (idx < f.dim0.size()) {
                return SurfaceSimplex{f.dim0[idx].v};
            }
            idx -= f.dim0.size();
            if (idx < f.dim1.size()) {
                return SurfaceSimplex{f.dim1[idx].v[0], f.dim1[idx].v[1]};
            }
            idx -= f.dim1.size();
            return SurfaceSimplex{f.dim2[idx].v[0], f.dim2[idx].v[1], f.dim2[idx].v[2]};
        })
        .method("get_simplex_value", [](const InterfaceSurface& s, int i) {
            const auto& f = s.filtration;
            if (i < 0 || i >= static_cast<int>(f.size())) {
                throw std::out_of_range("Simplex index out of range");
            }
            size_t idx = static_cast<size_t>(i);
            if (idx < f.dim0.size()) return f.dim0[idx].value;
            idx -= f.dim0.size();
            if (idx < f.dim1.size()) return f.dim1[idx].value;
            idx -= f.dim1.size();
            return f.dim2[idx].value;
        })
        .method("get_all_simplex_vertices_flat", [](const InterfaceSurface& s) {
            // [k_0, v_0_0, ..., k_1, v_1_0, ...] per filtration entry,
            // in the same entry order as get_all_simplex_values. Read
            // straight off the flat buckets (dim0|dim1|dim2 = entry order).
            const auto& f = s.filtration;
            std::vector<int32_t> result;
            result.reserve(f.dim0.size() * 2 + f.dim1.size() * 3 + f.dim2.size() * 4);
            for (const auto& e : f.dim0) {
                result.push_back(1);
                result.push_back(e.v);
            }
            for (const auto& e : f.dim1) {
                result.push_back(2);
                result.push_back(e.v[0]);
                result.push_back(e.v[1]);
            }
            for (const auto& e : f.dim2) {
                result.push_back(3);
                result.push_back(e.v[0]);
                result.push_back(e.v[1]);
                result.push_back(e.v[2]);
            }
            return result;
        })
        .method("get_all_simplex_values", [](const InterfaceSurface& s) {
            const auto& f = s.filtration;
            std::vector<double> result;
            result.reserve(f.size());
            for (const auto& e : f.dim0) result.push_back(e.value);
            for (const auto& e : f.dim1) result.push_back(e.value);
            for (const auto& e : f.dim2) result.push_back(e.value);
            return result;
        })
        .method("num_generating_tetrahedra", [](const InterfaceSurface& s) {
            return static_cast<int>(s.generating_simplices.tetrahedra.size());
        })
        .method("get_all_generating_tetrahedra", [](const InterfaceSurface& s) {
            std::vector<int> result;
            result.reserve(s.generating_simplices.tetrahedra.size() * 4);
            for (const auto& tet : s.generating_simplices.tetrahedra) {
                for (int v : tet) result.push_back(v);
            }
            return result;
        })
        .method("num_generating_free_triangles", [](const InterfaceSurface& s) {
            return static_cast<int>(s.generating_simplices.free_triangles.size());
        })
        .method("get_all_generating_free_triangles", [](const InterfaceSurface& s) {
            std::vector<int> result;
            for (const auto& tri : s.generating_simplices.free_triangles) {
                for (int v : tri) result.push_back(v);
            }
            return result;
        })
        .method("num_generating_free_edges", [](const InterfaceSurface& s) {
            return static_cast<int>(s.generating_simplices.free_edges.size());
        })
        .method("get_all_generating_free_edges", [](const InterfaceSurface& s) {
            std::vector<int> result;
            for (const auto& edge : s.generating_simplices.free_edges) {
                for (int v : edge) result.push_back(v);
            }
            return result;
        })
        .method("get_vertex_atom_indices_flat", [](const InterfaceSurface& s) {
            // Return as flat array: [n_atoms_0, atom_0_0, ..., n_atoms_1, atom_1_0, ...]
            std::vector<int> result;
            for (const auto& atoms : s.vertex_atom_indices) {
                result.push_back(static_cast<int>(atoms.size()));
                for (int a : atoms) result.push_back(a);
            }
            return result;
        });

    // InterfaceGenerator
    mod.add_type<InterfaceGenerator>("InterfaceGenerator")
        .constructor<>()
        .method("compute_interface_surface", [](
            InterfaceGenerator& gen,
            jlcxx::ArrayRef<double> flat_points,
            int n_points,
            jlcxx::ArrayRef<int> color_labels_arr,
            jlcxx::ArrayRef<double> radii_arr,
            bool alpha,
            bool lower_star
        ) {
            Points points = flat_points_to_cpp(flat_points, n_points);
            ColorLabels color_labels = julia_array_to_vector(color_labels_arr);
            Radii radii = julia_array_to_vector(radii_arr);
            return gen.compute_interface_surface(points, color_labels, radii, alpha, lower_star);
        })
        .method("compute_interface_surface", [](
            InterfaceGenerator& gen,
            jlcxx::ArrayRef<double> flat_points,
            int n_points,
            jlcxx::ArrayRef<int> color_labels_arr,
            double radius,
            bool lower_star
        ) {
            Points points = flat_points_to_cpp(flat_points, n_points);
            ColorLabels color_labels = julia_array_to_vector(color_labels_arr);
            return gen.compute_interface_surface(points, color_labels, radius, lower_star);
        })
        .method("get_multicolored_tetrahedra", [](
            InterfaceGenerator& gen,
            jlcxx::ArrayRef<double> flat_points,
            int n_points,
            jlcxx::ArrayRef<int> color_labels_arr,
            jlcxx::ArrayRef<double> radii_arr,
            bool alpha
        ) {
            Points points = flat_points_to_cpp(flat_points, n_points);
            ColorLabels color_labels = julia_array_to_vector(color_labels_arr);
            Radii radii = julia_array_to_vector(radii_arr);

            auto tets = gen.collect_multicolored_simplices(points, color_labels, radii, alpha).tetrahedra;

            // Convert to a flat array where each 4 consecutive integers represent one tetrahedron
            std::vector<int> result;
            result.reserve(tets.size() * 4);
            for (const auto& tet : tets) {
                result.push_back(tet[0]);
                result.push_back(tet[1]);
                result.push_back(tet[2]);
                result.push_back(tet[3]);
            }
            return result;
        });
}
