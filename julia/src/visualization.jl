"""
Visualization module for DelaunayInterfaces using GLMakie.

Provides functions to visualize interface surfaces and point clouds.
"""

using GLMakie
using GLMakie.Makie.GeometryBasics
using GLMakie.Colors

# =============================================================================
# Color Constants
#
# Configuration colors use Dark2_3 for distinguishing molecules/colors
# Interface colors use viridis for filtration values
# =============================================================================

const DEFAULT_INTERFACE_COLORMAP = :viridis
const DEFAULT_NUM_COLORS = 4

# Configuration/point cloud colormap - categorical colors for distinct labels
const CONF_GRADIENT = cgrad(:Dark2_4, 4, categorical=true)
const CONF_COLORMAP = [CONF_GRADIENT[i] for i in 1:4]
const DEFAULT_POINT_CLOUD_COLORMAP = :Dark2_4

# =============================================================================
# Mesh Generation
# =============================================================================

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

    # Extract vertex filtration values for coloring
    vertex_colors = Float64[]
    for (simplex, val) in surface.filtration
        if length(simplex) == 1
            push!(vertex_colors, val)
        end
    end

    return mesh, vertex_colors
end

# =============================================================================
# Interface Visualization
# =============================================================================

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

    if isempty(mesh_colors)
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
    draw_interface!(scene::LScene, surface::InterfaceSurface, points, color_labels, radii; kwargs...)

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
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
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
    interface_figure(surface, points, color_labels, radii; kwargs...)

Create a figure displaying an interface surface with point cloud overlay options.
"""
function interface_figure(
    surface::InterfaceSurface,
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
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

    draw_interface!(scene, surface, points, color_labels, radii;
        show_wireframe=show_wireframe,
        show_barycenters=show_barycenters,
        show_multicolored_points=show_multicolored_points,
        show_free_simplices=show_free_simplices,
        colormap=colormap,
        point_colormap=point_colormap
    )

    return fig
end

# =============================================================================
# Point Cloud Visualization
# =============================================================================

"""
    draw_point_cloud!(scene::LScene, points, color_labels, radii; kwargs...)

Draw a point cloud as spheres.

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

# =============================================================================
# Combined Visualization
# =============================================================================

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

    # Point cloud on the left
    gl_points = GridLayout(fig[1, 1])
    scene_points = LScene(gl_points[1, 1]; show_axis=show_axis)
    draw_point_cloud!(scene_points, points, color_labels, radii; colormap=point_colormap)

    # Interface on the right
    gl_interface = GridLayout(fig[1, 2])
    scene_interface = LScene(gl_interface[1, 1]; show_axis=show_axis, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))
    draw_interface!(scene_interface, surface;
        show_wireframe=show_wireframe,
        colormap=interface_colormap
    )

    if show_free_simplices
        draw_free_simplices!(scene_interface, surface)
    end

    # Add titles
    Label(gl_points[1, 1, Top()], "Atom Centers"; fontsize=16)
    Label(gl_interface[1, 1, Top()], "Interface"; fontsize=16)

    return fig
end

# =============================================================================
# Filtration Visualization
# =============================================================================

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

    # Observable data
    current_mesh = @lift(all_meshes[$slider.value])
    current_colors = @lift(all_colors[$slider.value])

    # Use final colorrange for consistent coloring
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

# =============================================================================
# Sequence Visualization
# =============================================================================

"""
    sequence_figure(surfaces, points_seq, labels_seq, radii_seq; kwargs...)

Create an interactive figure showing a sequence of interface surfaces.

# Arguments
- `surfaces::Vector{InterfaceSurface}`: Sequence of interface surfaces
- `points_seq::Vector{Vector{Vector{Float64}}}`: Points for each frame
- `labels_seq::Vector{Vector{Int}}`: Color labels for each frame
- `radii_seq::Vector{Vector{Float64}}`: Radii for each frame

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
    labels_seq::Vector{Vector{Int}},
    radii_seq::Vector{Vector{Float64}};
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

    individual_ranges = [(minimum(c), maximum(c)) for c in all_colors if !isempty(c)]

    fig = Figure(; fontsize=12)
    slider = Slider(fig[2, 1]; range=1:n, startvalue=1)
    scene = LScene(fig[1, 1];
        show_axis=false,
        scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],)
    )

    # Observable data
    current_mesh = @lift(all_meshes[$slider.value])
    current_colors = @lift(all_colors[$slider.value])

    colorrange = if global_colorrange
        all_vals = vcat(all_colors...)
        isempty(all_vals) ? (0.0, 1.0) : (minimum(all_vals), maximum(all_vals))
    else
        @lift(isempty(individual_ranges) ? (0.0, 1.0) : individual_ranges[$slider.value])
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
        current_free_verts = @lift(all_free_data[$slider.value][2])
        linesegments!(scene, current_free_edges; color=:orange, linewidth=2)
        scatter!(scene, current_free_verts; color=:red, markersize=10)
    end

    if show_multicolored_points
        current_points = @lift([Point3f(p...) for p in points_seq[$slider.value]])
        current_point_colors = @lift(
            cgrad(point_colormap, DEFAULT_NUM_COLORS; categorical=true)[labels_seq[$slider.value]]
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

# =============================================================================
# Helper Functions
# =============================================================================

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
    colors = cgrad(colormap, DEFAULT_NUM_COLORS; categorical=true)[color_labels]

    scatter!(scene, point_positions;
        color=colors,
        strokecolor=:black,
        strokewidth=1
    )
end

"""
    draw_free_simplices!(scene, surface; edge_color=:orange, point_color=:red, linewidth=2, markersize=10)

Draw the free simplices of the interface surface — filtration edges not part of any
triangle (as line segments between barycenters) and filtration vertices not part of
any edge (as points at barycenter positions).
"""
function draw_free_simplices!(
    scene::LScene,
    surface::InterfaceSurface;
    edge_color=:orange,
    point_color=:red,
    linewidth::Real=2,
    markersize::Real=10
)
    free_edge_pts, free_vert_pts = compute_free_simplex_data(surface)

    if !isempty(free_edge_pts)
        linesegments!(scene, free_edge_pts; color=edge_color, linewidth=linewidth)
    end

    if !isempty(free_vert_pts)
        scatter!(scene, free_vert_pts; color=point_color, markersize=markersize)
    end
end

"""
    compute_free_simplex_data(surface::InterfaceSurface)

Find filtration simplices that are not faces of higher-dimensional simplices.

# Returns
- `Tuple{Vector{Point3f}, Vector{Point3f}}`: (free edge endpoints for linesegments, free vertex positions)
"""
function compute_free_simplex_data(surface::InterfaceSurface)
    barycenters = [Point3f(v...) for v in surface.vertices]

    # Collect all edges and triangles from the filtration
    edges = Set{Tuple{Int32,Int32}}()
    triangle_edges = Set{Tuple{Int32,Int32}}()
    edge_verts = Set{Int32}()

    for (simplex, _) in surface.filtration
        if length(simplex) == 2
            e = minmax(simplex[1], simplex[2])
            push!(edges, e)
        elseif length(simplex) == 3
            for i in 1:3, j in (i+1):3
                push!(triangle_edges, minmax(simplex[i], simplex[j]))
            end
        end
    end

    # Free edges: in filtration but not face of any triangle
    free_edges = setdiff(edges, triangle_edges)
    free_edge_pts = Point3f[]
    for (i, j) in free_edges
        push!(free_edge_pts, barycenters[i], barycenters[j])
        push!(edge_verts, i, j)
    end

    # Also collect vertices that are endpoints of triangle edges
    for (i, j) in triangle_edges
        push!(edge_verts, i, j)
    end

    # Free vertices: in filtration but not endpoint of any edge
    all_edge_verts = edge_verts
    free_vert_pts = Point3f[]
    for (simplex, _) in surface.filtration
        if length(simplex) == 1 && !(simplex[1] in all_edge_verts)
            push!(free_vert_pts, barycenters[simplex[1]])
        end
    end

    return free_edge_pts, free_vert_pts
end

# =============================================================================
# Tetrahedron Visualization (Dual-Panel)
# =============================================================================

"""
    tetrahedron_subdivision_figure(points, colors, surface; title="")

Create a dual-panel figure for tetrahedron visualization:
- Left (LScene): Viridis-colored interface with wireframe, points, and tetrahedron edges
- Right (Axis3): Monocolor interface with wireframe and barycenter dots

# Arguments
- `points::Vector{Vector{Float64}}`: Tetrahedron vertices (4 points)
- `colors::Vector{Int}`: Color labels for each vertex
- `surface::InterfaceSurface`: The computed interface surface

# Keyword Arguments
- `title::String`: Figure title (default: "")

# Returns
- `Figure`: The GLMakie figure
"""
function tetrahedron_subdivision_figure(
    points::Vector{Vector{Float64}},
    colors::Vector{Int},
    surface::InterfaceSurface;
    title::String=""
)
    fig = Figure()

    # Left: LScene with viridis distance-colored interface
    Label(fig[1, 1, Top()], "Tetrahedron and Interface"; fontsize=16)
    scene_left = LScene(fig[1, 1]; show_axis=false)

    vertices = surface.vertices
    triangles = [Int.(simplex) for (simplex, _) in surface.filtration if length(simplex) == 3]

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Get filtration values for vertices (0-simplices) to color by distance
        vertex_vals = Dict{Int, Float64}()
        for (simplex, val) in surface.filtration
            if length(simplex) == 1
                vertex_vals[simplex[1]] = val
            end
        end
        mesh_colors = [get(vertex_vals, i, 0.0) for i in 1:length(vertices)]

        mesh!(scene_left, mesh_obj;
            color=mesh_colors,
            colormap=:viridis,
            colorrange=(minimum(mesh_colors), maximum(mesh_colors)),
            shading=NoShading
        )
        wireframe!(scene_left, mesh_obj; color=:white, linewidth=1)
    end

    # Draw original points colored by label
    pts_mat = reduce(hcat, points)'
    point_colors = [CONF_COLORMAP[mod1(c, 4)] for c in colors]
    scatter!(scene_left, pts_mat[:, 1], pts_mat[:, 2], pts_mat[:, 3];
             color=point_colors, markersize=20)

    # Draw tetrahedron edges with gradient colors
    n = length(points)
    edge_pts = Point3f[]
    edge_colors = RGBA[]
    for i in 1:n, j in (i+1):n
        p1, p2 = points[i], points[j]
        push!(edge_pts, Point3f(p1...), Point3f(p2...))
        c1 = CONF_COLORMAP[mod1(colors[i], 4)]
        c2 = CONF_COLORMAP[mod1(colors[j], 4)]
        push!(edge_colors, RGBA(c1), RGBA(c2))
    end
    linesegments!(scene_left, edge_pts; color=edge_colors, linewidth=2)

    # Right: Axis3 with monocolor mesh, wireframe, and vertices
    ax2 = Axis3(fig[1, 2]; aspect=:data, title="Interface and Barycenters")

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Mesh with transparency and NO shading
        mesh!(ax2, mesh_obj; color=RGBAf(0.27, 0.51, 0.71, 0.7), shading=NoShading)

        # Wireframe
        wireframe!(ax2, mesh_obj; color=:black, linewidth=1)

        # Red barycenter dots
        bary_pts = reduce(hcat, vertices)'
        scatter!(ax2, bary_pts[:, 1], bary_pts[:, 2], bary_pts[:, 3];
                 color=:red, markersize=12)
    end

    if !isempty(title)
        Label(fig[0, :], title; fontsize=20)
    end

    return fig
end

"""
    pointcloud_subdivision_figure(points, colors, surface; title="")

Create a dual-panel figure for point cloud visualization:
- Left (LScene): Viridis-colored interface with wireframe and points
- Right (Axis3): Monocolor interface with wireframe and barycenter dots (with lighting)

# Arguments
- `points::Vector{Vector{Float64}}`: Point cloud coordinates
- `colors::Vector{Int}`: Color labels for each point
- `surface::InterfaceSurface`: The computed interface surface

# Keyword Arguments
- `title::String`: Figure title (default: "")

# Returns
- `Figure`: The GLMakie figure
"""
function pointcloud_subdivision_figure(
    points::Vector{Vector{Float64}},
    colors::Vector{Int},
    surface::InterfaceSurface;
    title::String=""
)
    fig = Figure()

    # Left: LScene with viridis distance-colored interface + points
    Label(fig[1, 1, Top()], "Pointcloud and Interface"; fontsize=16)
    scene_left = LScene(fig[1, 1]; show_axis=false)

    vertices = surface.vertices
    triangles = [Int.(simplex) for (simplex, _) in surface.filtration if length(simplex) == 3]

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Get filtration values for vertices (0-simplices) to color by distance
        vertex_vals = Dict{Int, Float64}()
        for (simplex, val) in surface.filtration
            if length(simplex) == 1
                vertex_vals[simplex[1]] = val
            end
        end
        mesh_colors = [get(vertex_vals, i, 0.0) for i in 1:length(vertices)]

        mesh!(scene_left, mesh_obj;
            color=mesh_colors,
            colormap=:viridis,
            colorrange=(minimum(mesh_colors), maximum(mesh_colors)),
            shading=NoShading
        )
        wireframe!(scene_left, mesh_obj; color=:white, linewidth=1)
    end

    # Draw original points colored by label
    pts_mat = reduce(hcat, points)'
    point_colors = [CONF_COLORMAP[mod1(c, 4)] for c in colors]
    scatter!(scene_left, pts_mat[:, 1], pts_mat[:, 2], pts_mat[:, 3];
             color=point_colors, markersize=10)

    # Right: Axis3 with monocolor mesh, wireframe, and vertices (WITH lighting)
    ax2 = Axis3(fig[1, 2]; aspect=:data, title="Interface and Barycenters")

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Mesh with transparency and lighting enabled
        mesh!(ax2, mesh_obj; color=RGBAf(0.27, 0.51, 0.71, 0.7))

        # Wireframe
        wireframe!(ax2, mesh_obj; color=:black, linewidth=1)

        # Red barycenter dots
        bary_pts = reduce(hcat, vertices)'
        scatter!(ax2, bary_pts[:, 1], bary_pts[:, 2], bary_pts[:, 3];
                 color=:red, markersize=8)
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

    # LScene with lighting for the interface
    scene = LScene(fig[1, 1]; show_axis=false)

    vertices = surface.vertices
    triangles = [Int.(simplex) for (simplex, _) in surface.filtration if length(simplex) == 3]

    if !isempty(triangles) && !isempty(vertices)
        points_gb = [Point3f(v...) for v in vertices]
        faces_gb = [TriangleFace(t...) for t in triangles]
        mesh_obj = GeometryBasics.Mesh(points_gb, faces_gb)

        # Get filtration values for vertices to color by distance
        vertex_vals = Dict{Int, Float64}()
        for (simplex, val) in surface.filtration
            if length(simplex) == 1
                vertex_vals[simplex[1]] = val
            end
        end
        mesh_colors = [get(vertex_vals, i, 0.0) for i in 1:length(vertices)]

        # Mesh with viridis coloring (lighting enabled by default in LScene)
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

# =============================================================================
# Exports
# =============================================================================

export generate_colored_mesh
export draw_interface!, interface_figure
export draw_point_cloud!, point_cloud_figure
export draw_free_simplices!, compute_free_simplex_data
export interface_and_point_cloud_figure
export filtration_figure
export sequence_figure
export tetrahedron_subdivision_figure, pointcloud_subdivision_figure, interface_only_figure
export CONF_GRADIENT, CONF_COLORMAP
