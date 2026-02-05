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

# Load the Julia bindings
if !isdefined(Main, :DelaunayInterfaces)
    include(joinpath(@__DIR__, "..", "julia", "src", "DelaunayInterfaces.jl"))
end
using .DelaunayInterfaces

# =============================================================================
# Color Constants
# =============================================================================

const CONF_COLORS = [
    colorant"#1b9e77",  # teal
    colorant"#d95f02",  # orange
    colorant"#7570b3",  # purple
    colorant"#e7298a",  # pink
]

const INTERFACE_COLOR = RGBf(0.7, 0.7, 0.7)

# =============================================================================
# Tetrahedron Configurations
# =============================================================================

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

# =============================================================================
# Drawing Helper
# =============================================================================

function draw_tetrahedron_panel!(scene, points, colors, surface)
    # Draw order: tet edges -> mesh -> tet vertices -> interface edges

    # 1. Draw tetrahedron edges with gradient colors
    n = length(points)
    edge_points = Point3f[]
    edge_colors = RGBAf[]
    for i in 1:n, j in (i+1):n
        p1, p2 = Point3f(points[i]...), Point3f(points[j]...)
        c1 = CONF_COLORS[mod1(colors[i], 4)]
        c2 = CONF_COLORS[mod1(colors[j], 4)]
        push!(edge_points, p1, p2)
        push!(edge_colors, RGBAf(c1), RGBAf(c2))
    end
    linesegments!(scene, edge_points; color=edge_colors, linewidth=2)

    # 2. Draw interface mesh
    vertices = surface.vertices
    triangles = [Int.(simplex) for (simplex, _) in surface.filtration if length(simplex) == 3]

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Push mesh slightly back so lines render on top (depth_shift is in clip space)
        mesh!(scene, mesh_obj; color=INTERFACE_COLOR, shading=NoShading, depth_shift=0.0005f0)
    end

    # 3. Draw tetrahedron vertices as 3D spheres
    point_positions = [Point3f(p...) for p in points]
    point_colors = [CONF_COLORS[mod1(c, 4)] for c in colors]
    meshscatter!(scene, point_positions; color=point_colors, markersize=0.05, shading=NoShading)

    # 4. Draw interface triangle edges (last)
    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        tri_edge_points = Point3f[]
        for tri in triangles
            p1, p2, p3 = points_gb[tri[1]], points_gb[tri[2]], points_gb[tri[3]]
            push!(tri_edge_points, p1, p2, p2, p3, p3, p1)
        end
        linesegments!(scene, tri_edge_points; color=:black, linewidth=1.5)
    end
end

# =============================================================================
# Generate Figure
# =============================================================================

fig = Figure(size=(1000, 1000), backgroundcolor=:white)

positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

for (idx, (example, pos)) in enumerate(zip(EXAMPLES, positions))
    surface = InterfaceSurface(example.points, example.colors)

    # Use LScene for proper 3D rendering with depth testing
    scene = LScene(fig[pos...]; show_axis=false)

    draw_tetrahedron_panel!(scene, example.points, example.colors, surface)

    # Add title above each panel
    Label(fig[pos[1], pos[2], Top()], example.name; fontsize=16, padding=(0, 0, 5, 0))
end

# =============================================================================
# Interactive Display and Save
# =============================================================================

output_dir = joinpath(@__DIR__, "..", "assets")
mkpath(output_dir)
output_path = joinpath(output_dir, "tetrahedra_figure.png")

println("Interactive window opened.")
println("  - Rotate each panel to desired orientation")
println("  - Press 's' to save the figure")
println("  - Close the window when done")

# Add keyboard listener for saving
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press && event.key == Keyboard.s
        save(output_path, fig; px_per_unit=2)
        println("Figure saved to: $output_path")
    end
end

display(fig)

# Keep script running until window is closed
wait(GLMakie.Screen(fig.scene))
