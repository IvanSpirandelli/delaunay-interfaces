#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "delaunay_interfaces/interface_generation.hpp"

namespace py = pybind11;
using namespace delaunay_interfaces;

// Helper function to convert Nx3 numpy array to Points
Points numpy_to_points(const Eigen::Ref<const Eigen::MatrixXd>& arr) {
    Points points;
    points.reserve(arr.rows());
    for (Eigen::Index i = 0; i < arr.rows(); ++i) {
        points.push_back(Point3D(arr(i, 0), arr(i, 1), arr(i, 2)));
    }
    return points;
}

// Helper function to convert Points to Nx3 numpy array
Eigen::MatrixXd points_to_numpy(const Points& points) {
    Eigen::MatrixXd arr(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        arr(i, 0) = points[i].x();
        arr(i, 1) = points[i].y();
        arr(i, 2) = points[i].z();
    }
    return arr;
}

// Wrapper for InterfaceSurface with numpy-friendly vertices
struct InterfaceSurfacePy {
    Eigen::MatrixXd vertices;
    Filtration filtration;
    bool weighted;
    bool alpha;

    InterfaceSurfacePy(const InterfaceSurface& surface)
        : vertices(points_to_numpy(surface.vertices))
        , filtration(surface.filtration)
        , weighted(surface.weighted)
        , alpha(surface.alpha) {}
};

// Wrapper class with numpy-friendly interface
class InterfaceGeneratorPy {
    InterfaceGenerator gen_;
public:
    InterfaceGeneratorPy() = default;

    InterfaceSurfacePy compute_interface_surface(
        const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
        const ColorLabels& color_labels,
        const Radii& radii = {},
        bool weighted = true,
        bool alpha = true
    ) {
        Points points = numpy_to_points(points_arr);
        auto surface = gen_.compute_interface_surface(points, color_labels, radii, weighted, alpha);
        return InterfaceSurfacePy(surface);
    }

    Tetrahedra get_multicolored_tetrahedra(
        const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
        const ColorLabels& color_labels,
        const Radii& radii = {},
        bool weighted = true,
        bool alpha = true
    ) {
        Points points = numpy_to_points(points_arr);
        return gen_.get_multicolored_tetrahedra(points, color_labels, radii, weighted, alpha);
    }
};

// Wrapper for convenience function
std::pair<Eigen::MatrixXd, Filtration> get_barycentric_subdivision_and_filtration_py(
    const Eigen::Ref<const Eigen::MatrixXd>& points_arr,
    const ColorLabels& color_labels,
    const Radii& radii = {},
    bool weighted = true,
    bool alpha = true
) {
    Points points = numpy_to_points(points_arr);
    auto [vertices, filtration] = get_barycentric_subdivision_and_filtration(points, color_labels, radii, weighted, alpha);
    return {points_to_numpy(vertices), filtration};
}

PYBIND11_MODULE(delaunay_interfaces, m) {
    m.doc() = "DelaunayInterfaces: Compute interface surfaces from multicolored point clouds";

    // Bind InterfaceSurface (Python wrapper)
    py::class_<InterfaceSurfacePy>(m, "InterfaceSurface")
        .def_readonly("vertices", &InterfaceSurfacePy::vertices,
            "Barycenter vertices as Nx3 numpy array")
        .def_readonly("filtration", &InterfaceSurfacePy::filtration,
            "List of (simplex, filtration_value) tuples")
        .def_readonly("weighted", &InterfaceSurfacePy::weighted,
            "Whether weighted Delaunay/alpha complex was used")
        .def_readonly("alpha", &InterfaceSurfacePy::alpha,
            "Whether alpha complex was used");

    // Bind InterfaceGenerator (Python wrapper)
    py::class_<InterfaceGeneratorPy>(m, "InterfaceGenerator")
        .def(py::init<>())
        .def("compute_interface_surface", &InterfaceGeneratorPy::compute_interface_surface,
            py::arg("points"),
            py::arg("color_labels"),
            py::arg("radii") = Radii{},
            py::arg("weighted") = true,
            py::arg("alpha") = true,
            "Compute the interface surface from colored points\n\n"
            "Parameters\n"
            "----------\n"
            "points : numpy.ndarray (Nx3)\n"
            "    The input point cloud as Nx3 array\n"
            "color_labels : list of int\n"
            "    Color label for each point\n"
            "radii : list of float, optional\n"
            "    Radius for each point (required if weighted=True)\n"
            "weighted : bool, default=True\n"
            "    Use weighted Delaunay/alpha complex\n"
            "alpha : bool, default=True\n"
            "    Use alpha complex (vs Delaunay complex)\n\n"
            "Returns\n"
            "-------\n"
            "InterfaceSurface\n"
            "    The computed interface surface")
        .def("get_multicolored_tetrahedra", &InterfaceGeneratorPy::get_multicolored_tetrahedra,
            py::arg("points"),
            py::arg("color_labels"),
            py::arg("radii") = Radii{},
            py::arg("weighted") = true,
            py::arg("alpha") = true,
            "Get all multicolored tetrahedra from the complex\n\n"
            "Parameters\n"
            "----------\n"
            "points : numpy.ndarray (Nx3)\n"
            "    The input point cloud as Nx3 array\n"
            "color_labels : list of int\n"
            "    Color label for each point\n"
            "radii : list of float, optional\n"
            "    Radius for each point (required if weighted=True)\n"
            "weighted : bool, default=True\n"
            "    Use weighted Delaunay/alpha complex\n"
            "alpha : bool, default=True\n"
            "    Use alpha complex (vs Delaunay complex)\n\n"
            "Returns\n"
            "-------\n"
            "list of arrays\n"
            "    List of tetrahedra (each is array of 4 vertex indices)");

    // Convenience function
    m.def("get_barycentric_subdivision_and_filtration",
        &get_barycentric_subdivision_and_filtration_py,
        py::arg("points"),
        py::arg("color_labels"),
        py::arg("radii") = Radii{},
        py::arg("weighted") = true,
        py::arg("alpha") = true,
        "Compute barycentric subdivision and filtration\n\n"
        "Parameters\n"
        "----------\n"
        "points : numpy.ndarray (Nx3)\n"
        "    The input point cloud as Nx3 array\n"
        "color_labels : list of int\n"
        "    Color label for each point\n"
        "radii : list of float, optional\n"
        "    Radius for each point (required if weighted=True)\n"
        "weighted : bool, default=True\n"
        "    Use weighted Delaunay/alpha complex\n"
        "alpha : bool, default=True\n"
        "    Use alpha complex (vs Delaunay complex)\n\n"
        "Returns\n"
        "-------\n"
        "tuple of (vertices, filtration)\n"
        "    vertices: Nx3 numpy array of barycenter points\n"
        "    filtration: list of (simplex, filtration_value) tuples");

    // Version info
    m.attr("__version__") = "0.1.0";
}
