#include "jlcxx/jlcxx.hpp"
#include "delaunay_interfaces/interface_generation.hpp"
#include <vector>
#include <string>
#include <tuple>

using namespace delaunay_interfaces;

// Mark types as non-mirrored
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
    mod.method("version", []() { return std::string("0.1.0"); });

    // InterfaceSurface type
    mod.add_type<InterfaceSurface>("InterfaceSurfaceCxx")
        .method("num_vertices", [](const InterfaceSurface& s) {
            return static_cast<int>(s.vertices.size());
        })
        .method("num_simplices", [](const InterfaceSurface& s) {
            return static_cast<int>(s.filtration.size());
        })
        .method("is_weighted", [](const InterfaceSurface& s) {
            return s.weighted;
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
            if (i < 0 || i >= static_cast<int>(s.filtration.size())) {
                throw std::out_of_range("Simplex index out of range");
            }
            return std::get<0>(s.filtration[i]);
        })
        .method("get_simplex_value", [](const InterfaceSurface& s, int i) {
            if (i < 0 || i >= static_cast<int>(s.filtration.size())) {
                throw std::out_of_range("Simplex index out of range");
            }
            return std::get<1>(s.filtration[i]);
        })
        .method("num_generating_tetrahedra", [](const InterfaceSurface& s) {
            return static_cast<int>(s.simplices.generating_tetrahedra.size());
        })
        .method("get_all_generating_tetrahedra", [](const InterfaceSurface& s) {
            std::vector<int> result;
            result.reserve(s.simplices.generating_tetrahedra.size() * 4);
            for (const auto& tet : s.simplices.generating_tetrahedra) {
                for (int v : tet) result.push_back(v);
            }
            return result;
        })
        .method("num_generating_free_triangles", [](const InterfaceSurface& s) {
            return static_cast<int>(s.simplices.generating_free_triangles.size());
        })
        .method("get_all_generating_free_triangles", [](const InterfaceSurface& s) {
            std::vector<int> result;
            for (const auto& tri : s.simplices.generating_free_triangles) {
                for (int v : tri) result.push_back(v);
            }
            return result;
        })
        .method("num_generating_free_edges", [](const InterfaceSurface& s) {
            return static_cast<int>(s.simplices.generating_free_edges.size());
        })
        .method("get_all_generating_free_edges", [](const InterfaceSurface& s) {
            std::vector<int> result;
            for (const auto& edge : s.simplices.generating_free_edges) {
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
            bool weighted,
            bool alpha,
            bool lower_star
        ) {
            Points points = flat_points_to_cpp(flat_points, n_points);
            ColorLabels color_labels = julia_array_to_vector(color_labels_arr);
            Radii radii = julia_array_to_vector(radii_arr);
            return gen.compute_interface_surface(points, color_labels, radii, weighted, alpha, lower_star);
        })
        .method("get_multicolored_tetrahedra", [](
            InterfaceGenerator& gen,
            jlcxx::ArrayRef<double> flat_points,
            int n_points,
            jlcxx::ArrayRef<int> color_labels_arr,
            jlcxx::ArrayRef<double> radii_arr,
            bool weighted,
            bool alpha
        ) {
            Points points = flat_points_to_cpp(flat_points, n_points);
            ColorLabels color_labels = julia_array_to_vector(color_labels_arr);
            Radii radii = julia_array_to_vector(radii_arr);

            auto tets = gen.get_multicolored_tetrahedra(points, color_labels, radii, weighted, alpha);

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
