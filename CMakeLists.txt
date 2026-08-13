cmake_minimum_required(VERSION 3.15)
project(DelaunayInterfaces VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_JULIA_BINDINGS "Build Julia bindings (requires CxxWrap's JlCxx)" ON)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(delaunay_interfaces
    src/interface_generation.cpp
    src/barycentric_subdivision.cpp
    src/simplified_subdivision.cpp
)

target_include_directories(delaunay_interfaces PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delaunay_interfaces PUBLIC
    CGAL::CGAL
    Eigen3::Eigen
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(BUILD_JULIA_BINDINGS)
    find_package(JlCxx QUIET)
    if(JlCxx_FOUND)
        add_subdirectory(julia)
    else()
        message(STATUS "JlCxx not found - skipping Julia bindings. "
            "Point CMAKE_PREFIX_PATH at CxxWrap.prefix_path() to enable them.")
    endif()
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

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

include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/DelaunayInterfacesConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DelaunayInterfacesConfig.cmake
    INSTALL_DESTINATION lib/cmake/DelaunayInterfaces
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DelaunayInterfacesConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DelaunayInterfacesConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DelaunayInterfacesConfigVersion.cmake
    DESTINATION lib/cmake/DelaunayInterfaces
)
