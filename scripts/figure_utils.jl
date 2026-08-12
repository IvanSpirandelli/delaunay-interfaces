# Shared helpers for the figure-generation scripts.

using GLMakie

if !isdefined(Main, :DelaunayInterfaces)
    include(joinpath(@__DIR__, "..", "julia", "src", "DelaunayInterfaces.jl"))
end
using .DelaunayInterfaces

# Mesh/coloring helpers (generate_colored_mesh, vertex_filtration_values, ...)
include(joinpath(@__DIR__, "..", "julia", "src", "visualization.jl"))

"""
    interactive_save(fig, output_path)

Display `fig` in an interactive window; pressing 's' saves it to `output_path`.
Blocks until the window is closed.
"""
function interactive_save(fig::Figure, output_path::AbstractString)
    mkpath(dirname(output_path))

    println("Interactive window opened.")
    println("  - Rotate each panel to desired orientation")
    println("  - Press 's' to save the figure")
    println("  - Close the window when done")

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.s
            save(output_path, fig; px_per_unit=2)
            println("Figure saved to: $output_path")
        end
    end

    display(fig)
    wait(GLMakie.Screen(fig.scene))
end
