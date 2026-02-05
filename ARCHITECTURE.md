# DelaunayInterfaces Architecture

## Overview

C++ library for computing interface surfaces from multicolored 3D point clouds using Delaunay/alpha complexes and barycentric subdivision. Provides bindings for Python (pybind11) and Julia (CxxWrap.jl).

## Directory Structure

```
DelaunayInterfaces/
├── include/delaunay_interfaces/   # Public C++ API headers
│   ├── types.hpp                   # Type definitions and aliases
│   ├── interface_generation.hpp    # Main algorithm interface
│   ├── barycentric_subdivision.hpp # Subdivision logic
│   └── chromatic_partitioning.hpp  # Partitioning utilities
│
├── src/                     # C++ implementation
│   ├── interface_generation.cpp
│   └── barycentric_subdivision.cpp
│
├── python/                  # Python bindings (pybind11)
│   ├── bindings.cpp
│   └── setup.py
│
├── julia/                   # Julia bindings (CxxWrap.jl)
│   ├── bindings.cpp
│   └── DelaunayInterfaces.jl
│
├── examples/                # Example programs
└── tests/                   # Test suite
```

## Component Dependencies

```
┌─────────────────────────────────┐
│         CGAL + Eigen3           │
│   (Geometry + Linear Algebra)   │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│       C++ Core Library          │
│  (interface_generation.cpp)     │
│  (barycentric_subdivision.cpp)  │
└───────────────┬─────────────────┘
                │
        ┌───────┴───────┐
        ▼               ▼
┌──────────────┐ ┌──────────────┐
│ Python       │ │ Julia        │
│ (pybind11)   │ │ (CxxWrap.jl) │
└──────────────┘ └──────────────┘
```

## Core Modules

### interface_generation
Main entry point for computing interface surfaces.

| Function | Purpose |
|----------|---------|
| `generate_interface` | Compute interface from point cloud |
| `regular_delaunay_interface` | Regular Delaunay triangulation |
| `weighted_delaunay_interface` | Weighted Delaunay triangulation |
| `weighted_alpha_interface` | Weighted alpha complexes |

### barycentric_subdivision
Handles subdivision of multicolored tetrahedra into interface triangles.

| Function | Purpose |
|----------|---------|
| `compute_barycenter` | Compute barycenter of points |
| `extend_scaffold_*` | Scaffold extension for partition types |
| `get_filtration_value` | Distance-based filtration calculation |

### chromatic_partitioning (header-only)
Groups vertices by color for chromatic analysis. Inline utilities in `chromatic_partitioning.hpp`.

| Pattern | Description |
|---------|-------------|
| 2-2 | Two colors, two vertices each |
| 3-1 | Three of one color, one of another |
| 2-1-1 | Two-one-one distribution |
| 1-1-1-1 | Four different colors |

## Type System (types.hpp)

```cpp
using Point3 = Eigen::Vector3d;
using Points = std::vector<Point3>;
using Colors = std::vector<int>;
using Weights = std::vector<double>;
using Triangle = std::array<Point3, 3>;

struct ComplexConfig {
    bool use_weighted;
    bool use_alpha;
    double alpha_value;
};

struct InterfaceSurface {
    std::vector<Triangle> triangles;
    std::vector<double> filtration_values;
};
```

## Algorithm Flow

```
Input: Points + Colors [+ Weights]
            │
            ▼
    Delaunay/Alpha Complex (CGAL)
            │
            ▼
    Extract Multicolored Tetrahedra
    (2+ different vertex colors)
            │
            ▼
    Chromatic Partitioning
    (group by color distribution)
            │
            ▼
    Barycentric Subdivision
    (generate interface triangles)
            │
            ▼
Output: InterfaceSurface (triangles + filtration)
```

## Invariants

- All computations use Eigen for linear algebra (no raw arrays)
- CGAL types are aliased locally in implementation files only
- Python/Julia bindings expose identical API
- Filtration values are always computed during subdivision

## Layer Boundaries

- **Public API**: `include/delaunay_interfaces/*.hpp`
- **Implementation**: `src/*.cpp` (CGAL details hidden)
- **Bindings**: Wrap public API only, no CGAL exposure

## Dependencies

| Package | Purpose | Required |
|---------|---------|----------|
| CGAL | Computational geometry | Yes |
| Eigen3 | Linear algebra | Yes |
| pybind11 | Python bindings | Optional |
| CxxWrap.jl | Julia bindings | Optional |
