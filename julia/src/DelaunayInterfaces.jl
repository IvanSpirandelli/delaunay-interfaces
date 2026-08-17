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

# Re-export C++ types (0-based indices; prefer the Julia API below)
export InterfaceGenerator, InterfaceSurfaceCxx
export compute_interface_surface, get_multicolored_tetrahedra
export num_vertices, num_simplices, is_alpha
export get_vertex, get_all_vertices, get_simplex_vertices, get_simplex_value
export get_all_simplex_vertices_flat, get_all_simplex_values
export num_generating_tetrahedra, get_all_generating_tetrahedra
export num_generating_free_triangles, get_all_generating_free_triangles
export num_generating_free_edges, get_all_generating_free_edges
export get_vertex_atom_indices_flat, is_lower_star

# Julia API (1-based indices)
export InterfaceSurface
export get_triangles, get_edges, get_vertex_simplices

"""
    InterfaceSurface

Julia wrapper for interface surface results with convenient accessors.

# Fields
- `vertices::Vector{Vector{Float64}}` - Barycenter coordinates
- `filtration::Vector{Tuple{Vector{Int32}, Float64}}` - Simplices with filtration values
- `generating_tetrahedra::Matrix{Int}` - Generating multicolored tetrahedra, one per row
- `generating_free_triangles::Vector{Vector{Int}}` - Generating free multicolored triangles
- `generating_free_edges::Vector{Vector{Int}}` - Generating free multicolored edges
- `vertex_atom_indices::Vector{Vector{Int}}` - For each vertex, the sorted input atom indices it is the barycenter of
- `alpha::Bool` - Whether alpha complex was used
- `lower_star::Bool` - Whether lower star filtration was used
"""
struct InterfaceSurface
    vertices::Vector{Vector{Float64}}
    filtration::Vector{Tuple{Vector{Int32}, Float64}}
    generating_tetrahedra::Matrix{Int}
    generating_free_triangles::Vector{Vector{Int}}
    generating_free_edges::Vector{Vector{Int}}
    vertex_atom_indices::Vector{Vector{Int}}
    alpha::Bool
    lower_star::Bool
end

function InterfaceSurface(cxx_surface::InterfaceSurfaceCxx)
    # Extract vertices via one bulk FFI call (flat [x1,y1,z1,x2,...])
    n_verts = num_vertices(cxx_surface)
    flat_verts = collect(get_all_vertices(cxx_surface))
    vertices = [flat_verts[3i-2:3i] for i in 1:n_verts]

    # Extract filtration from two bulk arrays (convert 0-based C++ indices
    # to 1-based Julia indices)
    n_simplices = num_simplices(cxx_surface)
    flat_simplices = collect(get_all_simplex_vertices_flat(cxx_surface))
    simplex_values = collect(get_all_simplex_values(cxx_surface))
    filtration = Vector{Tuple{Vector{Int32}, Float64}}(undef, n_simplices)
    idx = 1
    for i in 1:n_simplices
        k = flat_simplices[idx]
        idx += 1
        filtration[i] = (flat_simplices[idx:idx+k-1] .+ Int32(1), simplex_values[i])
        idx += k
    end

    # Extract multicolored simplices (convert 0-based to 1-based)
    n_tets = num_generating_tetrahedra(cxx_surface)
    if n_tets > 0
        flat_tets = collect(get_all_generating_tetrahedra(cxx_surface))
        generating_tetrahedra = reshape(flat_tets, 4, :)' .+ 1
    else
        generating_tetrahedra = Matrix{Int}(undef, 0, 4)
    end

    n_free_tris = num_generating_free_triangles(cxx_surface)
    if n_free_tris > 0
        flat_tris = collect(get_all_generating_free_triangles(cxx_surface))
        generating_free_triangles = [flat_tris[i:i+2] .+ 1 for i in 1:3:length(flat_tris)]
    else
        generating_free_triangles = Vector{Vector{Int}}()
    end

    n_free_edg = num_generating_free_edges(cxx_surface)
    if n_free_edg > 0
        flat_edges = collect(get_all_generating_free_edges(cxx_surface))
        generating_free_edges = [flat_edges[i:i+1] .+ 1 for i in 1:2:length(flat_edges)]
    else
        generating_free_edges = Vector{Vector{Int}}()
    end

    # Extract vertex-to-atom mapping (convert 0-based to 1-based)
    flat_vai = collect(get_vertex_atom_indices_flat(cxx_surface))
    vertex_atom_indices = Vector{Vector{Int}}()
    idx = 1
    while idx <= length(flat_vai)
        n_atoms = flat_vai[idx]
        idx += 1
        push!(vertex_atom_indices, flat_vai[idx:idx+n_atoms-1] .+ 1)
        idx += n_atoms
    end

    InterfaceSurface(vertices, filtration, generating_tetrahedra, generating_free_triangles, generating_free_edges,
                     vertex_atom_indices, is_alpha(cxx_surface), is_lower_star(cxx_surface))
end

"""
    InterfaceSurface(points, color_labels[, radii]; alpha=!isempty(radii), lower_star=false)

Compute the interface surface from a colored point cloud.

# Arguments
- `points::Vector{Vector{Float64}}`: Vector of 3D points
- `color_labels::Vector{Int}`: Color label for each point (at least 2 distinct colors)
- `radii::Vector{Float64}`: Radius for each point; non-empty radii select the weighted complex
- `alpha::Bool`: Use alpha complex filtering (default: `true` iff radii given)
- `lower_star::Bool`: Use lower star filtration (max) instead of upper star (min) (default: `false`)

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
    alpha::Bool=!isempty(radii),
    lower_star::Bool=false
)
    gen = InterfaceGenerator()
    # Convert to CxxWrap-compatible types
    n_points = length(points)
    flat_points = reduce(vcat, points)
    color_labels_i32 = Int32.(color_labels)
    cxx_surface = compute_interface_surface(gen, flat_points, n_points, color_labels_i32, radii, alpha, lower_star)
    return InterfaceSurface(cxx_surface)
end

"""
    InterfaceSurface(points, color_labels, radius::Real; lower_star=false)

Uniform-radius variant: builds the alpha complex with parameter `radius^2`
(routed through CGAL's unweighted alpha shape). Always an alpha complex —
a radius has no effect on the plain Delaunay complex, so omit it for that.
"""
function InterfaceSurface(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radius::Real;
    lower_star::Bool=false
)
    gen = InterfaceGenerator()
    n_points = length(points)
    flat_points = reduce(vcat, points)
    color_labels_i32 = Int32.(color_labels)
    cxx_surface = compute_interface_surface(gen, flat_points, n_points, color_labels_i32, Float64(radius), lower_star)
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
    get_vertex_simplices(surface::InterfaceSurface)

Extract only the vertices (0-simplices) from the filtration.

# Returns
- `Vector{Tuple{Vector{Int32}, Float64}}`: Vertices with their filtration values
"""
function get_vertex_simplices(surface::InterfaceSurface)
    return filter(s -> length(s[1]) == 1, surface.filtration)
end

"""
    get_multicolored_tetrahedra(points, color_labels[, radii]; alpha=!isempty(radii))

Get all multicolored tetrahedra from the Delaunay/alpha complex.

# Returns
- `Matrix{Int}`: Matrix where each row is a tetrahedron with 4 vertex indices
"""
function get_multicolored_tetrahedra(
    points::Vector{Vector{Float64}},
    color_labels::Vector{Int},
    radii::Vector{Float64}=Float64[];
    alpha::Bool=!isempty(radii)
)
    gen = InterfaceGenerator()
    n_points = length(points)
    flat_points = reduce(vcat, points)
    color_labels_i32 = Int32.(color_labels)
    flat_result = get_multicolored_tetrahedra(gen, flat_points, n_points, color_labels_i32, radii, alpha)

    if isempty(flat_result)
        return Matrix{Int}(undef, 0, 4)
    end
    # Convert 0-based C++ indices to 1-based Julia indices
    return reshape(collect(flat_result), 4, :)' .+ 1
end

# Visualization module is included separately when GLMakie is available
# Use: include(joinpath(pkgdir(DelaunayInterfaces), "src", "visualization.jl"))
# Or simply `using DelaunayInterfaces; using GLMakie` and then include

end # module
