# Tests assert with <cassert>; keep assertions active even in Release builds
# (Release adds -DNDEBUG, which would compile every assert away).
add_compile_options(-UNDEBUG)

add_executable(test_basic test_basic.cpp)
target_link_libraries(test_basic PRIVATE delaunay_interfaces)

add_test(NAME BasicTest COMMAND test_basic)

# Combinatorial property tests (no external ground truth needed)
add_executable(test_combinatorial test_combinatorial.cpp)
target_link_libraries(test_combinatorial PRIVATE delaunay_interfaces)

add_test(NAME CombinatorialTest COMMAND test_combinatorial)

# Custom barycenter verification test
add_executable(test_custom_barycenter test_custom_barycenter.cpp)
target_link_libraries(test_custom_barycenter PRIVATE delaunay_interfaces)

add_test(NAME CustomBarycenterTest COMMAND test_custom_barycenter)

# Alpha complex construction tests
add_executable(test_alpha_construction test_alpha_construction.cpp)
target_link_libraries(test_alpha_construction PRIVATE delaunay_interfaces)
target_compile_definitions(test_alpha_construction PRIVATE
    TEST_OUTPUT_DIR="${CMAKE_CURRENT_BINARY_DIR}/failed-random-generated-tests"
)

add_test(NAME AlphaConstructionTest COMMAND test_alpha_construction)

# Simplified subdivision tests
add_executable(test_simplified test_simplified.cpp)
target_link_libraries(test_simplified PRIVATE delaunay_interfaces)

add_test(NAME SimplifiedTest COMMAND test_simplified)

# Benchmark (not a test: prints timings, asserts nothing; run manually)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
if(BUILD_BENCHMARKS)
    add_executable(benchmark_speed benchmark_speed.cpp)
    target_link_libraries(benchmark_speed PRIVATE delaunay_interfaces)
endif()
