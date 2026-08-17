"""
Visualization module for DelaunayInterfaces using GLMakie.

Provides functions to visualize interface surfaces and point clouds.
"""

using GLMakie
using GLMakie.Makie.GeometryBasics
using GLMakie.Colors

const DEFAULT_INTERFACE_COLORMAP = :viridis
const DEFAULT_NUM_COLORS = 4

# Categorical colors for point/molecule labels; viridis is used for filtration values.
const POINT_LABEL_GRADIENT = cgrad(:Dark2_4, DEFAULT_NUM_COLORS, categorical=true)
const POINT_LABEL_COLORS = [POINT_LABEL_GRADIENT[i] for i in 1:DEFAULT_NUM_COLORS]
const DEFAULT_POINT_CLOUD_COLORMAP = :Dark2_4

"""
    vertex_filtration_values(surface::InterfaceSurface)

Filtration value of each barycenter vertex, keyed by vertex ID.
"""
function vertex_filtration_values(surface::InterfaceSurface)
    vertex_vals = Dict{Int, Float64}()
    for (simplex, val) in surface.filtration
        if length(simplex) == 1
            vertex_vals[simplex[1]] = val
        end
    end
    return vertex_vals
end

"""
    generate_colored_mesh(surface::InterfaceSurface; max_value=Inf)

Generate a mesh and vertex colors from an interface surface.

# Arguments
- `surface::InterfaceSurface`: The interface surface to visualize
- `max_value::Real`: Maximum filtration value to include (default: `Inf`)

# Returns
- `Tuple{GeometryBasics.Mesh, Vector{Float64}}`: Mesh and per-vertex colors
"""
function generate_colored_mesh(surface::InterfaceSurface; max_value::Real=Inf)
    faces = TriangleFace{Int}[]
    for (simplex, val) in surface.filtration
        if length(simplex) == 3 && val <= max_value
            push!(faces, TriangleFace(simplex[1], simplex[2], simplex[3]))
        end
    end

    points = [Point3f(v...) for v in surface.vertices]
    mesh = GeometryBasics.Mesh(points, faces)

    vertex_vals = vertex_filtration_values(surface)
    vertex_colors = [get(vertex_vals, i, 0.0) for i in 1:length(surface.vertices)]

    return mesh, vertex_colors
end

"""
    draw_interface!(scene::LScene, surface::InterfaceSurface; kwargs...)

Draw an interface surface to an existing LScene.

# Keyword Arguments
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `show_barycenters::Bool`: Show barycenter points (default: `false`)
- `colormap`: Colormap for the surface (default: `:viridis`)
"""
function draw_interface!(
    scene::LScene,
    surface::InterfaceSurface;
    show_wireframe::Bool=false,
    show_barycenters::Bool=false,
    colormap=DEFAULT_INTERFACE_COLORMAP
)
    mesh_geom, mesh_colors = generate_colored_mesh(surface)

    if isempty(mesh_colors) || isempty(GeometryBasics.faces(mesh_geom))
        return
    end

    colorrange = (minimum(mesh_colors), maximum(mesh_colors))

    mesh!(scene, mesh_geom;
        color=mesh_colors,
        colorrange=colorrange,
        colormap=colormap,
        shading=NoShading
    )

    if show_wireframe
        wireframe!(scene, mesh_geom; color=:white, linewidth=1)
    end

    if show_barycenters
        draw_barycenters!(scene, [Point3f(v...) for v in surface.vertices])
    end
end

"""
    draw_interface!(scene::LScene, surface::InterfaceSurface, points, color_labels; kwargs...)

Draw an interface surface with optional multicolored point overlay.

# Keyword Arguments
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `show_barycenters::Bool`: Show barycenter points (default: `false`)
- `show_multicolored_points::Bool`: Show input points (default: `false`)
- `show_free_simplices::Bool`: Show free edges/vertices of the subdivision (default: `false`)
- `colormap`: Colormap for the surface (default: `:viridis`)
- `point_colormap`: Colormap for points (default: `:Dark2_4`)
"""
function draw_interface!(
    scene::LScene,
    surface::InterfaceSurface,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int};
    show_wireframe::Bool=false,
    show_barycenters::Bool=false,
    show_multicolored_points::Bool=false,
    show_free_simplices::Bool=false,
    colormap=DEFAULT_INTERFACE_COLORMAP,
    point_colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    draw_interface!(scene, surface;
        show_wireframe=show_wireframe,
        show_barycenters=show_barycenters,
        colormap=colormap
    )

    if show_multicolored_points
        draw_multicolored_points!(scene, points, color_labels; colormap=point_colormap)
    end

    if show_free_simplices
        draw_free_simplices!(scene, surface)
    end
end

"""
    interface_figure(surface::InterfaceSurface; kwargs...)

Create a figure displaying an interface surface.

# Keyword Arguments
- `show_axis::Bool`: Show coordinate axes (default: `false`)
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `show_barycenters::Bool`: Show barycenter points (default: `false`)
- `colormap`: Colormap for the surface (default: `:viridis`)

# Returns
- `Figure`: The GLMakie figure
"""
function interface_figure(
    surface::InterfaceSurface;
    show_axis::Bool=false,
    show_wireframe::Bool=false,
    show_barycenters::Bool=false,
    colormap=DEFAULT_INTERFACE_COLORMAP
)
    fig = Figure()
    scene = LScene(fig[1, 1]; show_axis=show_axis)

    draw_interface!(scene, surface;
        show_wireframe=show_wireframe,
        show_barycenters=show_barycenters,
        colormap=colormap
    )

    return fig
end

"""
    interface_figure(surface, points, color_labels; kwargs...)

Create a figure displaying an interface surface with point cloud overlay options.
"""
function interface_figure(
    surface::InterfaceSurface,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int};
    show_axis::Bool=false,
    show_wireframe::Bool=false,
    show_barycenters::Bool=false,
    show_multicolored_points::Bool=false,
    show_free_simplices::Bool=false,
    colormap=DEFAULT_INTERFACE_COLORMAP,
    point_colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    fig = Figure()
    scene = LScene(fig[1, 1]; show_axis=show_axis)

    draw_interface!(scene, surface, points, color_labels;
        show_wireframe=show_wireframe,
        show_barycenters=show_barycenters,
        show_multicolored_points=show_multicolored_points,
        show_free_simplices=show_free_simplices,
        colormap=colormap,
        point_colormap=point_colormap
    )

    return fig
end

"""
    draw_point_cloud!(scene::LScene, points, color_labels, radii; kwargs...)

Draw a point cloud as spheres. Radii set the sphere sizes when given.

# Keyword Arguments
- `colormap`: Colormap for the points (default: `:Dark2_4`)
"""
function draw_point_cloud!(
    scene::LScene,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    point_positions = [Point3f(p...) for p in points]
    marker_sizes = isempty(radii) ? fill(0.1f0, length(points)) : Float32.(radii)

    meshscatter!(scene, point_positions;
        markersize=marker_sizes,
        color=color_labels,
        colormap=colormap
    )
end

"""
    point_cloud_figure(points, color_labels, radii; kwargs...)

Create a figure displaying a point cloud.

# Keyword Arguments
- `show_axis::Bool`: Show coordinate axes (default: `false`)
- `colormap`: Colormap for the points (default: `:Dark2_4`)

# Returns
- `Figure`: The GLMakie figure
"""
function point_cloud_figure(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    show_axis::Bool=false,
    colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    fig = Figure()
    scene = LScene(fig[1, 1]; show_axis=show_axis)

    draw_point_cloud!(scene, points, color_labels, radii; colormap=colormap)

    return fig
end

"""
    interface_and_point_cloud_figure(surface, points, color_labels, radii; kwargs...)

Create a side-by-side figure with point cloud and interface surface.

# Keyword Arguments
- `show_axis::Bool`: Show coordinate axes (default: `false`)
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `show_free_simplices::Bool`: Show free edges/vertices of the subdivision (default: `false`)
- `interface_colormap`: Colormap for the interface (default: `:viridis`)
- `point_colormap`: Colormap for the points (default: `:Dark2_4`)

# Returns
- `Figure`: The GLMakie figure
"""
function interface_and_point_cloud_figure(
    surface::InterfaceSurface,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    show_axis::Bool=false,
    show_wireframe::Bool=false,
    show_free_simplices::Bool=false,
    interface_colormap=DEFAULT_INTERFACE_COLORMAP,
    point_colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    fig = Figure()

    gl_points = GridLayout(fig[1, 1])
    scene_points = LScene(gl_points[1, 1]; show_axis=show_axis)
    draw_point_cloud!(scene_points, points, color_labels, radii; colormap=point_colormap)

    gl_interface = GridLayout(fig[1, 2])
    scene_interface = LScene(gl_interface[1, 1]; show_axis=show_axis)
    draw_interface!(scene_interface, surface;
        show_wireframe=show_wireframe,
        colormap=interface_colormap
    )

    if show_free_simplices
        draw_free_simplices!(scene_interface, surface)
    end

    Label(gl_points[1, 1, Top()], "Atom Centers"; fontsize=16)
    Label(gl_interface[1, 1, Top()], "Interface"; fontsize=16)

    return fig
end

"""
    filtration_figure(surface::InterfaceSurface; kwargs...)

Create an interactive figure showing the filtration levels.

# Keyword Arguments
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `colormap`: Colormap for the surface (default: `:viridis`)

# Returns
- `Figure`: The GLMakie figure with a slider to explore filtration levels
"""
function filtration_figure(
    surface::InterfaceSurface;
    show_wireframe::Bool=false,
    colormap=DEFAULT_INTERFACE_COLORMAP
)
    levels = sort!(unique([val for (_, val) in surface.filtration]))
    meshes_and_colors = [generate_colored_mesh(surface; max_value=lvl) for lvl in levels]
    all_meshes = first.(meshes_and_colors)
    all_colors = last.(meshes_and_colors)

    fig = Figure(; fontsize=12)
    slider = Slider(fig[2, 1]; range=1:length(levels), startvalue=length(levels))
    scene = LScene(fig[1, 1];
        show_axis=false,
        scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],)
    )

    current_mesh = @lift(all_meshes[$slider.value])
    current_colors = @lift(all_colors[$slider.value])

    # The final level's colorrange keeps colors consistent across slider positions.
    final_colors = last(all_colors)
    colorrange = isempty(final_colors) ? (0.0, 1.0) : (minimum(final_colors), maximum(final_colors))

    mesh!(scene, current_mesh;
        color=current_colors,
        colorrange=colorrange,
        colormap=colormap
    )

    if show_wireframe
        wireframe!(scene, current_mesh; color=:white, linewidth=1)
    end

    return fig
end

"""
    sequence_figure(surfaces, points_seq, color_labels_seq; kwargs...)

Create an interactive figure showing a sequence of interface surfaces.

# Arguments
- `surfaces::Vector{InterfaceSurface}`: Sequence of interface surfaces
- `points_seq::Vector{Vector{Vector{Float64}}}`: Points for each frame
- `color_labels_seq::Vector{Vector{Int}}`: Color labels for each frame

# Keyword Arguments
- `show_wireframe::Bool`: Show wireframe overlay (default: `false`)
- `show_multicolored_points::Bool`: Show input points (default: `false`)
- `show_free_simplices::Bool`: Show free edges/vertices of the subdivision (default: `false`)
- `global_colorrange::Bool`: Use global color range across all frames (default: `false`)
- `interface_colormap`: Colormap for the interface (default: `:viridis`)
- `point_colormap`: Colormap for the points (default: `:Dark2_4`)

# Returns
- `Figure`: The GLMakie figure with a slider to explore the sequence
"""
function sequence_figure(
    surfaces::Vector{InterfaceSurface},
    points_seq::Vector{Vector{Vector{Float64}}},
    color_labels_seq::Vector{Vector{Int}};
    show_wireframe::Bool=false,
    show_multicolored_points::Bool=false,
    show_free_simplices::Bool=false,
    global_colorrange::Bool=false,
    interface_colormap=DEFAULT_INTERFACE_COLORMAP,
    point_colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    n = length(surfaces)
    colored_meshes = [generate_colored_mesh(s) for s in surfaces]
    all_meshes = first.(colored_meshes)
    all_colors = last.(colored_meshes)

    # One range per frame so frames with no triangles fall back cleanly.
    individual_ranges = [isempty(c) ? (0.0, 1.0) : (minimum(c), maximum(c)) for c in all_colors]

    fig = Figure(; fontsize=12)
    slider = Slider(fig[2, 1]; range=1:n, startvalue=1)
    scene = LScene(fig[1, 1];
        show_axis=false,
        scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],)
    )

    current_mesh = @lift(all_meshes[$slider.value])
    current_colors = @lift(all_colors[$slider.value])

    colorrange = if global_colorrange
        all_vals = vcat(all_colors...)
        isempty(all_vals) ? (0.0, 1.0) : (minimum(all_vals), maximum(all_vals))
    else
        @lift(individual_ranges[$slider.value])
    end

    mesh!(scene, current_mesh;
        color=current_colors,
        colorrange=colorrange,
        colormap=interface_colormap
    )

    if show_wireframe
        wireframe!(scene, current_mesh; color=:white, linewidth=1)
    end

    if show_free_simplices
        all_free_data = [compute_free_simplex_data(s) for s in surfaces]
        current_free_edges = @lift(all_free_data[$slider.value][1])
        current_free_edge_colors = @lift(all_free_data[$slider.value][2])
        current_free_verts = @lift(all_free_data[$slider.value][3])
        current_free_vert_colors = @lift(all_free_data[$slider.value][4])
        linesegments!(scene, current_free_edges; color=current_free_edge_colors, colormap=interface_colormap, colorrange=colorrange, linewidth=2)
        scatter!(scene, current_free_verts; color=current_free_vert_colors, colormap=interface_colormap, colorrange=colorrange, markersize=10)
    end

    if show_multicolored_points
        current_points = @lift([Point3f(p...) for p in points_seq[$slider.value]])
        current_point_colors = @lift(
            cgrad(point_colormap, DEFAULT_NUM_COLORS; categorical=true)[color_labels_seq[$slider.value]]
        )
        scatter!(scene, current_points;
            color=current_point_colors,
            markersize=15,
            strokewidth=1,
            strokecolor=:black
        )
    end

    return fig
end

"""
    subdivision_figure(points, color_labels, surface; title="", show_free_simplices=false)

Create a dual-panel figure for any colored point cloud:
- Left (LScene): Viridis-colored interface with wireframe, input points, and multicolored edges
- Right (Axis3): Monocolor interface with wireframe and barycenter dots

Works for tetrahedra, bipyramids, and large point clouds alike — edges are drawn
only for bicolored edges present in the alpha/Delaunay complex (not all pairs).

# Arguments
- `points::Vector{Vector{Float64}}`: Input point coordinates
- `color_labels::Vector{Int}`: Color labels for each point
- `surface::InterfaceSurface`: The computed interface surface

# Keyword Arguments
- `title::String`: Figure title (default: "")
- `show_free_simplices::Bool`: Show free edges and isolated vertices of the subdivision (default: `false`)
- `show_monocolored_edges::Bool`: Draw monocolored edges of the input (all same-colored point pairs); intended for single-tetrahedron inputs (default: `false`)

# Returns
- `Figure`: The GLMakie figure
"""
function subdivision_figure(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    surface::InterfaceSurface;
    title::String="",
    show_free_simplices::Bool=false,
    show_monocolored_edges::Bool=false
)
    fig = Figure()

    Label(fig[1, 1, Top()], "Point Cloud and Interface"; fontsize=16)
    scene_left = LScene(fig[1, 1]; show_axis=false)

    mesh_obj, mesh_colors = generate_colored_mesh(surface)
    has_mesh = !isempty(mesh_colors) && !isempty(GeometryBasics.faces(mesh_obj))

    if has_mesh
        mesh!(scene_left, mesh_obj;
            color=mesh_colors,
            colormap=:viridis,
            colorrange=(minimum(mesh_colors), maximum(mesh_colors)),
            shading=NoShading
        )
        wireframe!(scene_left, mesh_obj; color=:white, linewidth=1)
    end

    draw_multicolored_edges!(scene_left, points, color_labels, surface)

    if show_monocolored_edges
        draw_monocolored_edges!(scene_left, points, color_labels)
    end

    if show_free_simplices
        draw_free_simplices!(scene_left, surface)
    end

    pts_mat = reduce(hcat, points)'
    point_colors = [POINT_LABEL_COLORS[mod1(c, DEFAULT_NUM_COLORS)] for c in color_labels]
    scatter!(scene_left, pts_mat[:, 1], pts_mat[:, 2], pts_mat[:, 3];
             color=point_colors, markersize=15)

    ax2 = Axis3(fig[1, 2]; aspect=:data, title="Interface and Barycenters")

    if has_mesh
        mesh!(ax2, mesh_obj; color=RGBAf(0.27, 0.51, 0.71, 0.7), shading=NoShading)
        wireframe!(ax2, mesh_obj; color=:black, linewidth=1)

        bary_pts = reduce(hcat, surface.vertices)'
        scatter!(ax2, bary_pts[:, 1], bary_pts[:, 2], bary_pts[:, 3];
                 color=:red, markersize=10)
    end

    if !isempty(title)
        Label(fig[0, :], title; fontsize=20)
    end

    return fig
end

"""
    interface_only_figure(surface; title="", show_wireframe=false)

Create a single-panel figure showing only the interface surface with lighting.

# Arguments
- `surface::InterfaceSurface`: The computed interface surface

# Keyword Arguments
- `title::String`: Figure title (default: "")
- `show_wireframe::Bool`: Show wireframe overlay (default: false)

# Returns
- `Figure`: The GLMakie figure
"""
function interface_only_figure(
    surface::InterfaceSurface;
    title::String="",
    show_wireframe::Bool=false
)
    fig = Figure()
    scene = LScene(fig[1, 1]; show_axis=false)

    mesh_obj, mesh_colors = generate_colored_mesh(surface)

    if !isempty(mesh_colors) && !isempty(GeometryBasics.faces(mesh_obj))
        # Lighting stays enabled here (no NoShading), unlike draw_interface!.
        mesh!(scene, mesh_obj;
            color=mesh_colors,
            colormap=:viridis,
            colorrange=(minimum(mesh_colors), maximum(mesh_colors))
        )

        if show_wireframe
            wireframe!(scene, mesh_obj; color=:white, linewidth=1)
        end
    end

    if !isempty(title)
        Label(fig[0, :], title; fontsize=20)
    end

    return fig
end

"""
    draw_barycenters!(scene, points; markersize=15)

Draw barycenter points with index labels.
"""
function draw_barycenters!(scene::LScene, points::Vector{Point3f}; markersize::Int=15)
    scatter!(scene, points; color=:black, markersize=markersize, overdraw=true)
    for (i, pt) in enumerate(points)
        text!(scene, pt;
            text=string(i),
            fontsize=15,
            color=:red,
            overdraw=true,
            align=(:center, :center)
        )
    end
end

"""
    draw_multicolored_points!(scene, points, color_labels; kwargs...)

Draw input points with colors based on their labels.
"""
function draw_multicolored_points!(
    scene::LScene,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int};
    colormap=DEFAULT_POINT_CLOUD_COLORMAP
)
    point_positions = [Point3f(p...) for p in points]
    colors = cgrad(colormap, DEFAULT_NUM_COLORS; categorical=true)[mod1.(color_labels, DEFAULT_NUM_COLORS)]

    scatter!(scene, point_positions;
        color=colors,
        strokecolor=:black,
        strokewidth=1
    )
end

"""
    draw_free_simplices!(scene, surface; colormap=:viridis, colorrange=nothing, linewidth=2, markersize=10)

Draw the free simplices of the interface surface — filtration edges not part of any
triangle (as line segments between barycenters) and filtration vertices not part of
any edge (as points at barycenter positions).

Colors are derived from vertex filtration values using the same colormap as the surface.
"""
function draw_free_simplices!(
    scene::LScene,
    surface::InterfaceSurface;
    colormap=DEFAULT_INTERFACE_COLORMAP,
    colorrange=nothing,
    linewidth::Real=2,
    markersize::Real=10
)
    free_edge_pts, free_edge_colors, free_vert_pts, free_vert_colors = compute_free_simplex_data(surface)

    if isnothing(colorrange)
        all_vals = collect(values(vertex_filtration_values(surface)))
        colorrange = isempty(all_vals) ? (0.0, 1.0) : (minimum(all_vals), maximum(all_vals))
    end

    if !isempty(free_edge_pts)
        linesegments!(scene, free_edge_pts; color=free_edge_colors, colormap=colormap, colorrange=colorrange, linewidth=linewidth)
    end

    if !isempty(free_vert_pts)
        scatter!(scene, free_vert_pts; color=free_vert_colors, colormap=colormap, colorrange=colorrange, markersize=markersize)
    end
end

"""
    compute_free_simplex_data(surface::InterfaceSurface)

Find filtration simplices that are not faces of higher-dimensional simplices.

# Returns
- `Tuple{Vector{Point3f}, Vector{Float64}, Vector{Point3f}, Vector{Float64}}`:
  (free edge endpoints, edge endpoint filtration values, free vertex positions, vertex filtration values)
"""
function compute_free_simplex_data(surface::InterfaceSurface)
    barycenters = [Point3f(v...) for v in surface.vertices]
    vertex_vals = vertex_filtration_values(surface)

    edges = Set{Tuple{Int32,Int32}}()
    triangle_edges = Set{Tuple{Int32,Int32}}()
    edge_verts = Set{Int32}()

    for (simplex, _) in surface.filtration
        if length(simplex) == 2
            push!(edges, minmax(simplex[1], simplex[2]))
        elseif length(simplex) == 3
            for i in 1:3, j in (i+1):3
                push!(triangle_edges, minmax(simplex[i], simplex[j]))
            end
        end
    end

    # Free edges: in the filtration but not a face of any triangle.
    free_edges = setdiff(edges, triangle_edges)
    free_edge_pts = Point3f[]
    free_edge_colors = Float64[]
    for (i, j) in free_edges
        push!(free_edge_pts, barycenters[i], barycenters[j])
        push!(free_edge_colors, get(vertex_vals, Int(i), 0.0), get(vertex_vals, Int(j), 0.0))
        push!(edge_verts, i, j)
    end

    for (i, j) in triangle_edges
        push!(edge_verts, i, j)
    end

    # Free vertices: in the filtration but not an endpoint of any edge.
    free_vert_pts = Point3f[]
    free_vert_colors = Float64[]
    for (simplex, val) in surface.filtration
        if length(simplex) == 1 && !(simplex[1] in edge_verts)
            push!(free_vert_pts, barycenters[simplex[1]])
            push!(free_vert_colors, val)
        end
    end

    return free_edge_pts, free_edge_colors, free_vert_pts, free_vert_colors
end

"""
    draw_multicolored_edges!(scene, points, color_labels, surface; linewidth=2)

Draw the multicolored edges of the generating complex between input points.
Collects unique bicolored edges from generating tetrahedra, free triangles,
and free edges, then draws them as gradient-colored line segments.
"""
function draw_multicolored_edges!(
    scene::LScene,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    surface::InterfaceSurface;
    linewidth::Real=2
)
    pts = [Point3f(p...) for p in points]
    mc_edges = Set{Tuple{Int,Int}}()

    for row in eachrow(surface.generating_tetrahedra)
        verts = collect(row)
        for i in 1:4, j in (i+1):4
            color_labels[verts[i]] != color_labels[verts[j]] && push!(mc_edges, minmax(verts[i], verts[j]))
        end
    end

    for tri in surface.generating_free_triangles
        for i in 1:3, j in (i+1):3
            color_labels[tri[i]] != color_labels[tri[j]] && push!(mc_edges, minmax(tri[i], tri[j]))
        end
    end

    for edge in surface.generating_free_edges
        push!(mc_edges, minmax(edge[1], edge[2]))
    end

    isempty(mc_edges) && return

    edge_pts = Point3f[]
    edge_cols = RGBA[]
    for (i, j) in mc_edges
        push!(edge_pts, pts[i], pts[j])
        push!(edge_cols, RGBA(POINT_LABEL_COLORS[mod1(color_labels[i], DEFAULT_NUM_COLORS)]),
                         RGBA(POINT_LABEL_COLORS[mod1(color_labels[j], DEFAULT_NUM_COLORS)]))
    end
    linesegments!(scene, edge_pts; color=edge_cols, linewidth=linewidth)
end

"""
    draw_monocolored_edges!(scene, points, color_labels; linewidth=2)

Draw the monocolored edges of the input complex: every pair of input points
sharing the same color label, drawn in that color. Intended for small inputs
(single tetrahedra, bipyramids) where all point pairs are edges of the input.
"""
function draw_monocolored_edges!(
    scene::LScene,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int};
    linewidth::Real=2
)
    pts = [Point3f(p...) for p in points]

    edge_pts = Point3f[]
    edge_cols = RGBA[]
    for i in 1:length(pts), j in (i+1):length(pts)
        color_labels[i] == color_labels[j] || continue
        push!(edge_pts, pts[i], pts[j])
        c = RGBA(POINT_LABEL_COLORS[mod1(color_labels[i], DEFAULT_NUM_COLORS)])
        push!(edge_cols, c, c)
    end

    isempty(edge_pts) && return
    linesegments!(scene, edge_pts; color=edge_cols, linewidth=linewidth)
end

export vertex_filtration_values, generate_colored_mesh
export draw_interface!, interface_figure
export draw_point_cloud!, point_cloud_figure
export draw_free_simplices!, compute_free_simplex_data
export draw_multicolored_edges!, draw_monocolored_edges!
export interface_and_point_cloud_figure
export filtration_figure
export sequence_figure
export subdivision_figure, interface_only_figure
export POINT_LABEL_GRADIENT, POINT_LABEL_COLORS
