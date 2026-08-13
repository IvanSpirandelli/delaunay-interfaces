#!/usr/bin/env julia
"""
delaunay-interfaces Julia Example

This example demonstrates how to use the Julia bindings for the
delaunay-interfaces C++ library.
"""

using DelaunayInterfaces

function main()
    println("delaunay-interfaces Julia Example")
    println("=" ^ 50)
    println()

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
    radii = fill(0.5, length(points))

    println("Input:")
    println("  Points: ", length(points))
    println("  Colors: 2 groups")
    println()

    # Example 1: Using InterfaceSurface constructor (unweighted Delaunay)
    println("Example 1: Unweighted Delaunay Interface")
    println("-" ^ 50)

    surface = InterfaceSurface(points, colors)

    println("  Barycenters: ", length(surface.vertices))
    println("  Filtration simplices: ", length(surface.filtration))
    println("  Weighted: ", surface.weighted)
    println("  Alpha: ", surface.alpha)
    println()

    # Show some filtration details
    if length(surface.filtration) > 0
        println("First 10 filtration simplices:")
        for (i, (simplex, value)) in enumerate(surface.filtration[1:min(10, end)])
            println("  ", i, ": dim=", length(simplex)-1, ", value=", round(value, digits=6))
        end
        println()
    end

    # Example 2: Get triangles and edges
    println("Example 2: Extract triangles and edges")
    println("-" ^ 50)

    triangles = get_triangles(surface)
    edges = get_edges(surface)
    vertex_simplices = get_vertices_simplices(surface)

    println("  Triangles (2-simplices): ", length(triangles))
    println("  Edges (1-simplices): ", length(edges))
    println("  Vertices (0-simplices): ", length(vertex_simplices))
    println()

    # Example 3: Compare different complex types
    println("Example 3: Compare different complex types")
    println("-" ^ 50)

    # Unweighted Delaunay
    surface_delaunay = InterfaceSurface(points, colors)
    n_del = length(surface_delaunay.vertices)
    t_del = length(get_triangles(surface_delaunay))
    println("  Unweighted Delaunay:  $n_del barycenters, $t_del triangles")

    # Weighted Delaunay
    surface_weighted = InterfaceSurface(points, colors, radii; weighted=true, alpha=false)
    n_wei = length(surface_weighted.vertices)
    t_wei = length(get_triangles(surface_weighted))
    println("  Weighted Delaunay:    $n_wei barycenters, $t_wei triangles")

    # Weighted Alpha
    surface_alpha = InterfaceSurface(points, colors, radii; weighted=true, alpha=true)
    n_alp = length(surface_alpha.vertices)
    t_alp = length(get_triangles(surface_alpha))
    println("  Weighted Alpha:       $n_alp barycenters, $t_alp triangles")

    println()

    # Example 4: Get multicolored tetrahedra
    println("Example 4: Get multicolored tetrahedra")
    println("-" ^ 50)

    mc_tets = get_multicolored_tetrahedra_wrapper(points, colors; weighted=false, alpha=false)

    println("  Number of multicolored tetrahedra: ", size(mc_tets, 1))
    if size(mc_tets, 1) > 0
        println("  First tetrahedron indices: ", mc_tets[1, :])
    end
    println()

    println("Done!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
