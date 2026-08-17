#!/usr/bin/env julia
"""Benchmark the Julia FFI boundary: C++ compute time vs. wrapper conversion.

Times the raw CxxWrap call (compute_interface_surface) separately from the
InterfaceSurface(cxx_surface) conversion, which crosses the boundary per
element. All-distinct colors maximize output size, stressing the conversion.

Usage: julia --project=julia scripts/benchmark_julia_boundary.jl [out.json] [sizes_csv] [repeats]
"""

using Random
using Statistics
using JSON
using DelaunayInterfaces

function bench(n::Int, repeats::Int)
    Random.seed!(1000003 * n)
    points = [rand(3) for _ in 1:n]
    flat_points = reduce(vcat, points)
    color_labels = Int32.(1:n)
    gen = InterfaceGenerator()

    cxx_ms = Float64[]
    wrap_ms = Float64[]
    n_vertices = 0
    n_simplices = 0
    for _ in 1:repeats
        t0 = time_ns()
        cxx_surface = compute_interface_surface(gen, flat_points, n, color_labels, Float64[], false, false)
        t1 = time_ns()
        surface = InterfaceSurface(cxx_surface)
        t2 = time_ns()
        push!(cxx_ms, (t1 - t0) / 1e6)
        push!(wrap_ms, (t2 - t1) / 1e6)
        n_vertices = length(surface.vertices)
        n_simplices = length(surface.filtration)
    end

    return Dict(
        "n" => n,
        "cxx_ms" => median(cxx_ms),
        "wrapper_ms" => median(wrap_ms),
        "n_vertices" => n_vertices,
        "n_simplices" => n_simplices,
    )
end

function main()
    out_path = length(ARGS) >= 1 ? ARGS[1] : ""
    sizes = length(ARGS) >= 2 ? parse.(Int, split(ARGS[2], ',')) : [1000, 10000, 30000]
    repeats = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 3

    bench(100, 1)  # warm-up: JIT-compile the timed path before measuring

    results = [bench(n, repeats) for n in sizes]

    println("Julia FFI boundary (median of $repeats repeats):")
    println("  N pts   |   cxx ms | wrapper ms")
    for r in results
        println("  ", lpad(r["n"], 7), " | ", lpad(round(r["cxx_ms"], digits=1), 8),
                " | ", lpad(round(r["wrapper_ms"], digits=1), 10))
    end

    if !isempty(out_path)
        open(out_path, "w") do io
            JSON.print(io, results, 1)
        end
        println("JSON written to $out_path")
    end
end

main()
