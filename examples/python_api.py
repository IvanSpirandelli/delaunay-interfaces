#!/usr/bin/env python3
"""
delaunay-interfaces Python Example

This example demonstrates how to use the Python bindings for the
delaunay-interfaces C++ library.
"""

import sys
import os

# Add build directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)
sys.path.insert(0, os.path.join(repo_root, 'build', 'python'))

try:
    import delaunay_interfaces as di
except ImportError:
    print("Error: delaunay_interfaces module not found.")
    print("Please build the Python bindings first:")
    print("  cd build && cmake .. && make")
    sys.exit(1)


def main():
    print("delaunay-interfaces Python Example")
    print("=" * 50)
    print()

    # Create a simple point cloud with two color groups
    points = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.5, 1.0],
        [2.0, 0.0, 0.0],
        [2.5, 1.0, 0.0],
        [2.5, 0.5, 1.0],
        [1.5, 0.5, 0.5]
    ]

    colors = [1, 1, 1, 1, 2, 2, 2, 2]
    radii = [0.5] * len(points)

    print("Input:")
    print(f"  Points: {len(points)}")
    print(f"  Color groups: {len(set(colors))}")
    print()

    gen = di.InterfaceGenerator()

    # Example 1: Unweighted Delaunay Interface
    print("Example 1: Unweighted Delaunay Interface")
    print("-" * 50)

    surface = gen.compute_interface_surface(points, colors)

    print(f"  Barycenters: {len(surface.vertices)}")
    print(f"  Filtration simplices: {len(surface.filtration)}")
    print(f"  Alpha: {surface.alpha}")
    print()

    # Show some filtration details
    if len(surface.filtration) > 0:
        print("First 10 filtration simplices:")
        for i, (simplex, value) in enumerate(surface.filtration[:10]):
            print(f"  {i+1}: dim={len(simplex)-1}, value={value:.6f}")
        print()

    # Example 2: Extract triangles and edges
    print("Example 2: Extract triangles and edges")
    print("-" * 50)

    triangles = [s for s, v in surface.filtration if len(s) == 3]
    edges = [s for s, v in surface.filtration if len(s) == 2]
    vertices = [s for s, v in surface.filtration if len(s) == 1]

    print(f"  Triangles (2-simplices): {len(triangles)}")
    print(f"  Edges (1-simplices): {len(edges)}")
    print(f"  Vertices (0-simplices): {len(vertices)}")
    print()

    # Example 3: Compare different complex types
    print("Example 3: Compare different complex types")
    print("-" * 50)

    # Unweighted Delaunay (no radii)
    surface_del = gen.compute_interface_surface(points, colors)
    n_del = len(surface_del.vertices)
    t_del = sum(1 for s, v in surface_del.filtration if len(s) == 3)
    print(f"  Unweighted Delaunay:  {n_del} barycenters, {t_del} triangles")

    # Weighted Delaunay
    surface_wei = gen.compute_interface_surface(points, colors, radii, alpha=False)
    n_wei = len(surface_wei.vertices)
    t_wei = sum(1 for s, v in surface_wei.filtration if len(s) == 3)
    print(f"  Weighted Delaunay:    {n_wei} barycenters, {t_wei} triangles")

    # Weighted Alpha (default when radii are given)
    surface_alp = gen.compute_interface_surface(points, colors, radii)
    n_alp = len(surface_alp.vertices)
    t_alp = sum(1 for s, v in surface_alp.filtration if len(s) == 3)
    print(f"  Weighted Alpha:       {n_alp} barycenters, {t_alp} triangles")

    # Uniform Alpha: a single radius builds the alpha complex with parameter
    # radius^2 (same surface as uniform per-point radii, cheaper to construct)
    surface_uni = gen.compute_interface_surface(points, colors, radius=0.5)
    n_uni = len(surface_uni.vertices)
    t_uni = sum(1 for s, v in surface_uni.filtration if len(s) == 3)
    print(f"  Uniform Alpha:        {n_uni} barycenters, {t_uni} triangles")

    print()
    print("Done!")


if __name__ == "__main__":
    main()
