find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Try to find pybind11
find_package(pybind11 CONFIG)

if(NOT pybind11_FOUND)
    # If not found, fetch it from GitHub
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(delaunay_interfaces_py
    bindings.cpp
)

target_link_libraries(delaunay_interfaces_py PRIVATE
    delaunay_interfaces
)

# Set output name to match Python import
set_target_properties(delaunay_interfaces_py PROPERTIES
    OUTPUT_NAME delaunay_interfaces
)

# Installation
install(TARGETS delaunay_interfaces_py
    LIBRARY DESTINATION ${Python3_SITEARCH}
)
