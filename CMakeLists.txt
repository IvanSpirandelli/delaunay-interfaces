cmake_minimum_required(VERSION 3.15)
project(DelaunayInterfaces VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_JULIA_BINDINGS "Build Julia bindings" ON)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Find required packages
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

# Fetch nlohmann/json for testing
include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

# Core library
add_library(delaunay_interfaces SHARED
    src/interface_generation.cpp
    src/barycentric_subdivision.cpp
)

target_include_directories(delaunay_interfaces PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delaunay_interfaces PUBLIC
    CGAL::CGAL
    Eigen3::Eigen
)

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

# Julia bindings
if(BUILD_JULIA_BINDINGS)
    add_subdirectory(julia)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Installation
install(TARGETS delaunay_interfaces
    EXPORT DelaunayInterfacesTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/delaunay_interfaces
    DESTINATION include
)

install(EXPORT DelaunayInterfacesTargets
    FILE DelaunayInterfacesTargets.cmake
    NAMESPACE DelaunayInterfaces::
    DESTINATION lib/cmake/DelaunayInterfaces
)
