#!/usr/bin/env julia
"""Generate README figure: protein dimer with point cloud and interface surface.

Usage:
  1. Run this script to open an interactive GLMakie window
  2. Rotate each panel to the desired orientation
  3. Press 's' to save the figure, or close the window
"""

using JSON
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
]

# =============================================================================
# Load Protein Data
# =============================================================================

data_path = joinpath(@__DIR__, "..", "tests", "data", "ground_truth_4bmg_dimer_alpha.json")
protein_data = JSON.parsefile(data_path)

inp = protein_data["input"]
points = [Vector{Float64}(p) for p in inp["points"]]
colors = Vector{Int}(inp["color_labels"])
radii = Vector{Float64}(inp["radii"])

println("Protein: 4bmg_dimer")
println("  Points: $(length(points))")
println("  Unique colors: $(length(unique(colors)))")

# =============================================================================
# Compute Interface
# =============================================================================

surface = InterfaceSurface(points, colors, radii; weighted=true, alpha=true)

num_triangles = count(x -> length(x[1]) == 3, surface.filtration)
println("Interface surface:")
println("  Vertices: $(length(surface.vertices))")
println("  Triangles: $num_triangles")

# =============================================================================
# Create Figure
# =============================================================================

fig = Figure(size=(1400, 700), backgroundcolor=:white)

# Probe radius used in surface generation
rs = 1.4
display_radii = radii .- rs

# Left panel: Point cloud (atoms as spheres)
scene_left = LScene(fig[1, 1]; show_axis=false)

point_positions = [Point3f(p...) for p in points]
point_colors = [CONF_COLORS[mod1(c, 2)] for c in colors]
marker_sizes = Float32.(display_radii)

meshscatter!(scene_left, point_positions;
    markersize=marker_sizes,
    color=point_colors
)

Label(fig[1, 1, Top()], "Atom Centers"; fontsize=18, padding=(0, 0, 10, 0))

# Right panel: Interface surface
scene_right = LScene(fig[1, 2]; show_axis=false)

vertices = surface.vertices
triangles = [Int.(simplex) for (simplex, _) in surface.filtration if length(simplex) == 3]

if !isempty(triangles) && !isempty(vertices)
    points_gb = [Point3f(v...) for v in vertices]
    faces_gb = [TriangleFace(t...) for t in triangles]
    mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

    # Get filtration values for coloring
    vertex_vals = Dict{Int, Float64}()
    for (simplex, val) in surface.filtration
        if length(simplex) == 1
            vertex_vals[simplex[1]] = val
        end
    end
    mesh_colors = [get(vertex_vals, i, 0.0) for i in 1:length(vertices)]

    mesh!(scene_right, mesh_obj;
        color=mesh_colors,
        colormap=:viridis,
        colorrange=(minimum(mesh_colors), maximum(mesh_colors)),
        shading=NoShading
    )
end

Label(fig[1, 2, Top()], "Interface Surface"; fontsize=18, padding=(0, 0, 10, 0))

# =============================================================================
# Interactive Display and Save
# =============================================================================

output_dir = joinpath(@__DIR__, "..", "assets")
mkpath(output_dir)
output_path = joinpath(output_dir, "protein_figure.png")

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
