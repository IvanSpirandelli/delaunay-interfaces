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

include(joinpath(@__DIR__, "figure_utils.jl"))

const CONF_COLORS = [
    colorant"#1b9e77",  # teal
    colorant"#d95f02",  # orange
]

data_path = joinpath(@__DIR__, "..", "tests", "data", "ground_truth_4bmg_dimer_delaunay.json")
protein_data = JSON.parsefile(data_path)

inp = protein_data["input"]
points = [Vector{Float64}(p) for p in inp["points"]]
colors = Vector{Int}(inp["color_labels"])
radii = Vector{Float64}(inp["radii"])

println("Protein: 4bmg_dimer")
println("  Points: $(length(points))")
println("  Unique colors: $(length(unique(colors)))")

surface = InterfaceSurface(points, colors, radii; weighted=true, alpha=true)

num_triangles = count(x -> length(x[1]) == 3, surface.filtration)
println("Interface surface:")
println("  Vertices: $(length(surface.vertices))")
println("  Triangles: $num_triangles")

fig = Figure(size=(1400, 700), backgroundcolor=:white)

# The radii include the probe (water) radius; subtract it back out for display.
rs = 1.4
display_radii = radii .- rs

scene_left = LScene(fig[1, 1]; show_axis=false)

meshscatter!(scene_left, [Point3f(p...) for p in points];
    markersize=Float32.(display_radii),
    color=[CONF_COLORS[mod1(c, 2)] for c in colors]
)

Label(fig[1, 1, Top()], "Atom Centers"; fontsize=18, padding=(0, 0, 10, 0))

scene_right = LScene(fig[1, 2]; show_axis=false)

mesh_obj, mesh_colors = generate_colored_mesh(surface)
if !isempty(mesh_colors) && !isempty(GeometryBasics.faces(mesh_obj))
    mesh!(scene_right, mesh_obj;
        color=mesh_colors,
        colormap=:viridis,
        colorrange=(minimum(mesh_colors), maximum(mesh_colors)),
        shading=NoShading
    )
end

Label(fig[1, 2, Top()], "Interface Surface"; fontsize=18, padding=(0, 0, 10, 0))

interactive_save(fig, joinpath(@__DIR__, "..", "assets", "protein_figure.png"))
