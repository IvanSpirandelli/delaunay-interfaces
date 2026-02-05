module DelaunayInterfaces

using CxxWrap

# Load the C++ library - use absolute path for reliable loading during precompilation
const LIB_PATH = let
    # @__DIR__ is julia/src, go up to project root then into build/julia
    project_root = dirname(dirname(@__DIR__))
    joinpath(project_root, "build", "julia", "libdelaunay_interfaces_jl")
end

@wrapmodule(() -> LIB_PATH)

function __init__()
    @initcxx
end

# Re-export C++ types
export InterfaceGenerator, InterfaceSurfaceCxx
export compute_interface_surface, get_multicolored_tetrahedra
export num_vertices, num_simplices, is_weighted, is_alpha
export get_vertex, get_all_vertices, get_simplex_vertices, get_simplex_value

"""
    InterfaceSurface

Julia wrapper for interface surface results with convenient accessors.

# Fields (accessible via methods)
- `vertices::Vector{Vector{Float64}}` - Barycenter coordinates
- `filtration::Vector{Tuple{Vector{Int32}, Float64}}` - Simplices with filtration values
- `weighted::Bool` - Whether weighted complex was used
- `alpha::Bool` - Whether alpha complex was used
"""
struct InterfaceSurface
    _cxx::InterfaceSurfaceCxx
    vertices::Vector{Vector{Float64}}
    filtration::Vector{Tuple{Vector{Int32}, Float64}}
    weighted::Bool
    alpha::Bool
end

function InterfaceSurface(cxx_surface::InterfaceSurfaceCxx)
    # Extract vertices
    n_verts = num_vertices(cxx_surface)
    vertices = Vector{Vector{Float64}}(undef, n_verts)
    for i in 0:(n_verts - 1)
        vertices[i + 1] = collect(get_vertex(cxx_surface, i))
    end

    # Extract filtration (convert 0-based C++ indices to 1-based Julia indices)
    n_simplices = num_simplices(cxx_surface)
    filtration = Vector{Tuple{Vector{Int32}, Float64}}(undef, n_simplices)
    for i in 0:(n_simplices - 1)
        simplex_verts = collect(get_simplex_vertices(cxx_surface, i)) .+ Int32(1)
        simplex_val = get_simplex_value(cxx_surface, i)
        filtration[i + 1] = (simplex_verts, simplex_val)
    end

    InterfaceSurface(cxx_surface, vertices, filtration, is_weighted(cxx_surface), is_alpha(cxx_surface))
end

"""
    InterfaceSurface(points, color_labels[, radii]; weighted=true, alpha=true)

Compute the interface surface from a colored point cloud.

# Arguments
- `points::Vector{Vector{Float64}}`: Vector of 3D points
- `color_labels::Vector{Int}`: Color label for each point (at least 2 distinct colors)
- `radii::Vector{Float64}`: Radius for each point (required if `weighted=true`)
- `weighted::Bool`: Use weighted Delaunay/alpha complex (default: `true`)
- `alpha::Bool`: Use alpha complex filtering (default: `true`)

# Returns
- `InterfaceSurface`: Object containing vertices and filtration data

# Examples
```julia
points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
colors = [1, 1, 2, 2]
radii = [0.5, 0.5, 0.5, 0.5]

surface = InterfaceSurface(points, colors, radii)
println("Vertices: ", length(surface.vertices))
println("Simplices: ", length(surface.filtration))
```
"""
function InterfaceSurface(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    weighted::Bool=!isempty(radii),  # Default: weighted if radii provided
    alpha::Bool=!isempty(radii)       # Default: alpha if radii provided
)
    gen = InterfaceGenerator()
    # Convert to CxxWrap-compatible types
    n_points = length(points)
    flat_points = reduce(vcat, points)
    color_labels_i32 = Int32.(color_labels)
    cxx_surface = compute_interface_surface(gen, flat_points, n_points, color_labels_i32, radii, weighted, alpha)
    return InterfaceSurface(cxx_surface)
end

"""
    get_triangles(surface::InterfaceSurface)

Extract only the triangle faces (2-simplices) from the filtration.

# Returns
- `Vector{Tuple{Vector{Int32}, Float64}}`: Triangles with their filtration values
"""
function get_triangles(surface::InterfaceSurface)
    return filter(s -> length(s[1]) == 3, surface.filtration)
end

"""
    get_edges(surface::InterfaceSurface)

Extract only the edges (1-simplices) from the filtration.

# Returns
- `Vector{Tuple{Vector{Int32}, Float64}}`: Edges with their filtration values
"""
function get_edges(surface::InterfaceSurface)
    return filter(s -> length(s[1]) == 2, surface.filtration)
end

"""
    get_vertices_simplices(surface::InterfaceSurface)

Extract only the vertices (0-simplices) from the filtration.

# Returns
- `Vector{Tuple{Vector{Int32}, Float64}}`: Vertices with their filtration values
"""
function get_vertices_simplices(surface::InterfaceSurface)
    return filter(s -> length(s[1]) == 1, surface.filtration)
end

"""
    get_multicolored_tetrahedra_wrapper(points, color_labels[, radii]; weighted=true, alpha=true)

Get all multicolored tetrahedra from the Delaunay/alpha complex.

# Returns
- `Matrix{Int}`: Matrix where each row is a tetrahedron with 4 vertex indices
"""
function get_multicolored_tetrahedra_wrapper(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    weighted::Bool=true,
    alpha::Bool=true
)
    gen = InterfaceGenerator()
    # Convert to CxxWrap-compatible types
    n_points = length(points)
    flat_points = reduce(vcat, points)
    color_labels_i32 = Int32.(color_labels)
    flat_result = get_multicolored_tetrahedra(gen, flat_points, n_points, color_labels_i32, radii, weighted, alpha)

    if isempty(flat_result)
        return Matrix{Int}(undef, 0, 4)
    end
    # Convert 0-based C++ indices to 1-based Julia indices
    return reshape(collect(flat_result), 4, :)' .+ 1
end

export InterfaceSurface
export get_triangles, get_edges, get_vertices_simplices
export get_multicolored_tetrahedra_wrapper

# Visualization module is included separately when GLMakie is available
# Use: include(joinpath(pkgdir(DelaunayInterfaces), "src", "visualization.jl"))
# Or simply `using DelaunayInterfaces; using GLMakie` and then include

end # module
