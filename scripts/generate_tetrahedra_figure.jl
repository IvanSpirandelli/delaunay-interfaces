#!/usr/bin/env julia
"""Generate README figure: 2x2 grid of tetrahedron partition examples.

Usage:
  1. Run this script to open an interactive GLMakie window
  2. Rotate each panel to the desired orientation
  3. Press 's' to save the figure, or close the window
"""

using LinearAlgebra
using GLMakie
using GeometryBasics

include(joinpath(@__DIR__, "figure_utils.jl"))

const INTERFACE_COLOR = RGBf(0.7, 0.7, 0.7)

function rotate_points(points, axis, angle)
    axis = normalize(axis)
    K = [0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0]
    R = I + sin(angle) * K + (1 - cos(angle)) * K^2
    [Vector{Float64}(R * p) for p in points]
end

function center_points(points)
    centroid = sum(points) / length(points)
    [p - centroid for p in points]
end

const UNIT_TETRAHEDRON_BASE = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
]
const UNIT_TETRAHEDRON = center_points(rotate_points(UNIT_TETRAHEDRON_BASE, [1.0, 0.5, 0.7], π/6))

const SHARED_TETRAHEDRON = [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, -1.0], [0.0, 1.0, 1.0]]

# 2-2: rotate around x-axis to show the interface plane better
const SHARED_22 = center_points(rotate_points(SHARED_TETRAHEDRON, [1.0, 0.0, 0.0], π/4))

# 2-1-1: tilt toward screen (rotate around x-axis)
const SHARED_211 = center_points(rotate_points(SHARED_TETRAHEDRON, [1.0, 0.0, 0.0], -π/6))

# 3-1: centered
const TET_31 = let r = 1.0
    center_points([[r, -1.0, 0.0], [-r/2, -1.0, r*sqrt(3)/2], [-r/2, -1.0, -r*sqrt(3)/2], [0.0, 1.0, 0.0]])
end

const EXAMPLES = [
    (name="3-1 Partition", points=TET_31, colors=[1, 1, 1, 2]),
    (name="2-2 Partition", points=SHARED_22, colors=[1, 1, 2, 2]),
    (name="2-1-1 Partition", points=SHARED_211, colors=[1, 1, 2, 3]),
    (name="1-1-1-1 Partition", points=UNIT_TETRAHEDRON, colors=[1, 2, 3, 4]),
]

function draw_tetrahedron_panel!(scene, points, colors, surface)
    # Draw order matters: tet edges, then mesh, then tet vertices, then interface edges.

    n = length(points)
    edge_points = Point3f[]
    edge_colors = RGBAf[]
    for i in 1:n, j in (i+1):n
        push!(edge_points, Point3f(points[i]...), Point3f(points[j]...))
        push!(edge_colors, RGBAf(POINT_LABEL_COLORS[mod1(colors[i], 4)]),
                           RGBAf(POINT_LABEL_COLORS[mod1(colors[j], 4)]))
    end
    linesegments!(scene, edge_points; color=edge_colors, linewidth=2)

    mesh_obj, mesh_colors = generate_colored_mesh(surface)
    has_mesh = !isempty(mesh_colors) && !isempty(GeometryBasics.faces(mesh_obj))

    if has_mesh
        # Push mesh slightly back so lines render on top (depth_shift is in clip space)
        mesh!(scene, mesh_obj; color=INTERFACE_COLOR, shading=NoShading, depth_shift=0.0005f0)
    end

    meshscatter!(scene, [Point3f(p...) for p in points];
        color=[POINT_LABEL_COLORS[mod1(c, 4)] for c in colors],
        markersize=0.05, shading=NoShading)

    if has_mesh
        verts = GeometryBasics.coordinates(mesh_obj)
        tri_edge_points = Point3f[]
        for tri in GeometryBasics.faces(mesh_obj)
            p1, p2, p3 = verts[tri[1]], verts[tri[2]], verts[tri[3]]
            push!(tri_edge_points, p1, p2, p2, p3, p3, p1)
        end
        linesegments!(scene, tri_edge_points; color=:black, linewidth=1.5)
    end
end

fig = Figure(size=(1000, 1000), backgroundcolor=:white)

positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

for (example, pos) in zip(EXAMPLES, positions)
    surface = InterfaceSurface(example.points, example.colors)

    scene = LScene(fig[pos...]; show_axis=false)
    draw_tetrahedron_panel!(scene, example.points, example.colors, surface)
    Label(fig[pos[1], pos[2], Top()], example.name; fontsize=16, padding=(0, 0, 5, 0))
end

interactive_save(fig, joinpath(@__DIR__, "..", "assets", "tetrahedra_figure.png"))
