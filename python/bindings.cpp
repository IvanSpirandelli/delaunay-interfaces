#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <optional>
#include "delaunay_interfaces/interface_generation.hpp"

namespace py = pybind11;
using namespace delaunay_interfaces;

Points numpy_to_points(const Eigen::Ref<const Eigen::MatrixXd>& arr) {
    Points points;
    points.reserve(arr.rows());
    for (Eigen::Index i = 0; i < arr.rows(); ++i) {
        points.push_back(Point3D(arr(i, 0), arr(i, 1), arr(i, 2)));
    }
    return points;
}

Eigen::MatrixXd points_to_numpy(const Points& points) {
    Eigen::MatrixXd arr(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        arr(i, 0) = points[i].x();
        arr(i, 1) = points[i].y();
        arr(i, 2) = points[i].z();
    }
    return arr;
}

struct InterfaceSurfacePy {
    Eigen::MatrixXd vertices;
    Filtration filtration;
    Tetrahedra generating_tetrahedra;
    std::vector<FreeTriangle> generating_free_triangles;
    std::vector<FreeEdge> generating_free_edges;
    VertexAtomIndices vertex_atom_indices;
    bool weighted;
    bool alpha;
    bool lower_star;

    InterfaceSurfacePy(const InterfaceSurface& surface)
        : vertices(points_to_numpy(surface.vertices))
        , filtration(surface.filtration)
        , generating_tetrahedra(surface.generating_simplices.generating_tetrahedra)
        , generating_free_triangles(surface.generating_simplices.generating_free_triangles)
        , generating_free_edges(surface.generating_simplices.generating_free_edges)
        , vertex_atom_indices(surface.vertex_atom_indices)
        , weighted(surface.weighted)
        , alpha(surface.alpha)
        , lower_star(surface.lower_star) {}
};

class InterfaceGeneratorPy {
    InterfaceGenerator gen_;
public:
    InterfaceGeneratorPy() = default;

    InterfaceSurfacePy compute_interface_surface(
        const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
        const ColorLabels& color_labels,
        const Radii& radii,
        std::optional<bool> weighted,
        std::optional<bool> alpha,
        bool lower_star
    ) {
        Points points = numpy_to_points(points_arr);
        auto surface = gen_.compute_interface_surface(
            points, color_labels, radii,
            weighted.value_or(!radii.empty()), alpha.value_or(!radii.empty()), lower_star);
        return InterfaceSurfacePy(surface);
    }

    Tetrahedra get_multicolored_tetrahedra(
        const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
        const ColorLabels& color_labels,
        const Radii& radii,
        std::optional<bool> weighted,
        std::optional<bool> alpha
    ) {
        Points points = numpy_to_points(points_arr);
        return gen_.collect_multicolored_simplices(
            points, color_labels, radii,
            weighted.value_or(!radii.empty()), alpha.value_or(!radii.empty())).generating_tetrahedra;
    }
};

struct SimplifiedSurfacePy {
    Eigen::MatrixXd vertices;
    std::vector<std::array<int32_t, 3>> triangles;
    std::vector<std::array<int32_t, 4>> quads;
    std::vector<std::array<int32_t, 2>> edges;
    VertexAtomIndices vertex_atom_indices;
    std::vector<double> vertex_filtration;
    Tetrahedra generating_tetrahedra;
    std::vector<FreeTriangle> generating_free_triangles;
    std::vector<FreeEdge> generating_free_edges;
    bool weighted;
    bool alpha;

    SimplifiedSurfacePy(const SimplifiedSurface& surface)
        : vertices(points_to_numpy(surface.vertices))
        , triangles(surface.triangles)
        , quads(surface.quads)
        , edges(surface.edges)
        , vertex_atom_indices(surface.vertex_atom_indices)
        , vertex_filtration(surface.vertex_filtration)
        , generating_tetrahedra(surface.generating_simplices.generating_tetrahedra)
        , generating_free_triangles(surface.generating_simplices.generating_free_triangles)
        , generating_free_edges(surface.generating_simplices.generating_free_edges)
        , weighted(surface.weighted)
        , alpha(surface.alpha) {}
};

SimplifiedSurfacePy compute_simplified_surface_py(
    const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> weighted,
    std::optional<bool> alpha
) {
    Points points = numpy_to_points(points_arr);
    auto surface = compute_simplified_surface(
        points, color_labels, radii,
        weighted.value_or(!radii.empty()), alpha.value_or(!radii.empty()));
    return SimplifiedSurfacePy(surface);
}

std::pair<Eigen::MatrixXd, Filtration> compute_barycentric_subdivision_and_filtration_py(
    const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
    const ColorLabels& color_labels,
    const Radii& radii,
    std::optional<bool> weighted,
    std::optional<bool> alpha,
    bool lower_star
) {
    Points points = numpy_to_points(points_arr);
    auto [vertices, filtration, simplices, vertex_atom_indices] = compute_barycentric_subdivision_and_filtration(
        points, color_labels, radii,
        weighted.value_or(!radii.empty()), alpha.value_or(!radii.empty()), lower_star);
    return {points_to_numpy(vertices), filtration};
}

PYBIND11_MODULE(delaunay_interfaces, m) {
    m.doc() = "DelaunayInterfaces: Compute interface surfaces from multicolored point clouds";

    py::class_<InterfaceSurfacePy>(m, "InterfaceSurface")
        .def_readonly("vertices", &InterfaceSurfacePy::vertices,
            "Barycenter vertices as Nx3 numpy array")
        .def_readonly("filtration", &InterfaceSurfacePy::filtration,
            "List of (simplex, filtration_value) tuples")
        .def_readonly("generating_tetrahedra", &InterfaceSurfacePy::generating_tetrahedra,
            "Generating multicolored tetrahedra (list of 4-element arrays)")
        .def_readonly("generating_free_triangles", &InterfaceSurfacePy::generating_free_triangles,
            "Generating free multicolored triangles not covered by tetrahedra")
        .def_readonly("generating_free_edges", &InterfaceSurfacePy::generating_free_edges,
            "Generating free multicolored edges not covered by tetrahedra or triangles")
        .def_readonly("vertex_atom_indices", &InterfaceSurfacePy::vertex_atom_indices,
            "For each barycenter vertex, the sorted list of input atom indices it is the barycenter of")
        .def_readonly("weighted", &InterfaceSurfacePy::weighted,
            "Whether weighted Delaunay/alpha complex was used")
        .def_readonly("alpha", &InterfaceSurfacePy::alpha,
            "Whether alpha complex was used")
        .def_readonly("lower_star", &InterfaceSurfacePy::lower_star,
            "Whether lower star filtration was used (max instead of min)");

    py::class_<InterfaceGeneratorPy>(m, "InterfaceGenerator")
        .def(py::init<>())
        .def("compute_interface_surface", &InterfaceGeneratorPy::compute_interface_surface,
            py::arg("points"),
            py::arg("color_labels"),
            py::arg("radii") = Radii{},
            py::arg("weighted") = py::none(),
            py::arg("alpha") = py::none(),
            py::arg("lower_star") = false,
            "Compute the interface surface from colored points\n\n"
            "Parameters\n"
            "----------\n"
            "points : numpy.ndarray (Nx3)\n"
            "    The input point cloud as Nx3 array\n"
            "color_labels : list of int\n"
            "    Color label for each point\n"
            "radii : list of float, optional\n"
            "    Radius for each point (required if weighted=True)\n"
            "weighted : bool, optional\n"
            "    Use weighted Delaunay/alpha complex (default: True iff radii given)\n"
            "alpha : bool, optional\n"
            "    Use alpha complex vs Delaunay complex (default: True iff radii given)\n"
            "lower_star : bool, default=False\n"
            "    Use lower star filtration (max) instead of upper star (min)\n\n"
            "Returns\n"
            "-------\n"
            "InterfaceSurface\n"
            "    The computed interface surface")
        .def("get_multicolored_tetrahedra", &InterfaceGeneratorPy::get_multicolored_tetrahedra,
            py::arg("points"),
            py::arg("color_labels"),
            py::arg("radii") = Radii{},
            py::arg("weighted") = py::none(),
            py::arg("alpha") = py::none(),
            "Get all multicolored tetrahedra from the complex\n\n"
            "Parameters\n"
            "----------\n"
            "points : numpy.ndarray (Nx3)\n"
            "    The input point cloud as Nx3 array\n"
            "color_labels : list of int\n"
            "    Color label for each point\n"
            "radii : list of float, optional\n"
            "    Radius for each point (required if weighted=True)\n"
            "weighted : bool, optional\n"
            "    Use weighted Delaunay/alpha complex (default: True iff radii given)\n"
            "alpha : bool, optional\n"
            "    Use alpha complex vs Delaunay complex (default: True iff radii given)\n\n"
            "Returns\n"
            "-------\n"
            "list of arrays\n"
            "    List of tetrahedra (each is array of 4 vertex indices)");

    m.def("compute_barycentric_subdivision_and_filtration",
        &compute_barycentric_subdivision_and_filtration_py,
        py::arg("points"),
        py::arg("color_labels"),
        py::arg("radii") = Radii{},
        py::arg("weighted") = py::none(),
        py::arg("alpha") = py::none(),
        py::arg("lower_star") = false,
        "Compute barycentric subdivision and filtration\n\n"
        "Parameters\n"
        "----------\n"
        "points : numpy.ndarray (Nx3)\n"
        "    The input point cloud as Nx3 array\n"
        "color_labels : list of int\n"
        "    Color label for each point\n"
        "radii : list of float, optional\n"
        "    Radius for each point (required if weighted=True)\n"
        "weighted : bool, optional\n"
        "    Use weighted Delaunay/alpha complex (default: True iff radii given)\n"
        "alpha : bool, optional\n"
        "    Use alpha complex vs Delaunay complex (default: True iff radii given)\n"
        "lower_star : bool, default=False\n"
        "    Use lower star filtration (max) instead of upper star (min)\n\n"
        "Returns\n"
        "-------\n"
        "tuple of (vertices, filtration)\n"
        "    vertices: Nx3 numpy array of barycenter points\n"
        "    filtration: list of (simplex, filtration_value) tuples");

    py::class_<SimplifiedSurfacePy>(m, "SimplifiedSurface")
        .def_readonly("vertices", &SimplifiedSurfacePy::vertices,
            "Atom-pair midpoint vertices as Nx3 numpy array")
        .def_readonly("triangles", &SimplifiedSurfacePy::triangles,
            "Triangle faces from 3-1 tets (list of 3-element arrays)")
        .def_readonly("quads", &SimplifiedSurfacePy::quads,
            "Quad faces from 2-2 tets (list of 4-element arrays, cyclic order)")
        .def_readonly("edges", &SimplifiedSurfacePy::edges,
            "Free edges from 2-1 free triangles (list of 2-element arrays)")
        .def_readonly("vertex_atom_indices", &SimplifiedSurfacePy::vertex_atom_indices,
            "For each vertex, the pair [atom_a, atom_b] of input atom indices")
        .def_readonly("vertex_filtration", &SimplifiedSurfacePy::vertex_filtration,
            "Filtration value (inter-atom distance) per vertex")
        .def_readonly("generating_tetrahedra", &SimplifiedSurfacePy::generating_tetrahedra,
            "Generating multicolored tetrahedra")
        .def_readonly("generating_free_triangles", &SimplifiedSurfacePy::generating_free_triangles,
            "Generating free multicolored triangles")
        .def_readonly("generating_free_edges", &SimplifiedSurfacePy::generating_free_edges,
            "Generating free multicolored edges")
        .def_readonly("weighted", &SimplifiedSurfacePy::weighted)
        .def_readonly("alpha", &SimplifiedSurfacePy::alpha);

    m.def("compute_simplified_surface",
        &compute_simplified_surface_py,
        py::arg("points"),
        py::arg("color_labels"),
        py::arg("radii") = Radii{},
        py::arg("weighted") = py::none(),
        py::arg("alpha") = py::none(),
        "Compute simplified 2-color interface surface\n\n"
        "Each vertex = midpoint of one cross-color atom pair.\n"
        "Requires exactly 2 distinct colors in the input.\n\n"
        "Parameters\n"
        "----------\n"
        "points : numpy.ndarray (Nx3)\n"
        "color_labels : list of int (exactly 2 distinct values)\n"
        "radii : list of float, optional\n"
        "weighted : bool, optional (default: True iff radii given)\n"
        "alpha : bool, optional (default: True iff radii given)\n\n"
        "Returns\n"
        "-------\n"
        "SimplifiedSurface");

    m.attr("__version__") = "0.1.0";
}
