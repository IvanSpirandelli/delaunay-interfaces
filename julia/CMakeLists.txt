find_package(JlCxx REQUIRED)

add_library(delaunay_interfaces_jl SHARED
    bindings.cpp
)

target_link_libraries(delaunay_interfaces_jl
    PUBLIC
        delaunay_interfaces
        JlCxx::cxxwrap_julia
)

# Installation
install(TARGETS delaunay_interfaces_jl
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
