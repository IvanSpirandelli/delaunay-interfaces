# Basic test support
enable_testing()

# Add a simple test executable
add_executable(test_basic test_basic.cpp)
target_link_libraries(test_basic PRIVATE delaunay_interfaces)

add_test(NAME BasicTest COMMAND test_basic)

# Ground truth validation test (small tetrahedra)
add_executable(test_ground_truth test_ground_truth.cpp)
target_link_libraries(test_ground_truth PRIVATE delaunay_interfaces)

add_test(NAME GroundTruthTest COMMAND test_ground_truth)

# Large examples ground truth test (proteins and random point clouds)
add_executable(test_large_examples test_large_examples.cpp)
target_link_libraries(test_large_examples PRIVATE delaunay_interfaces nlohmann_json::nlohmann_json)
target_compile_definitions(test_large_examples PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)

add_test(NAME LargeExamplesTest COMMAND test_large_examples)

# Combinatorial property tests (no external ground truth needed)
add_executable(test_combinatorial test_combinatorial.cpp)
target_link_libraries(test_combinatorial PRIVATE delaunay_interfaces)

add_test(NAME CombinatorialTest COMMAND test_combinatorial)

# Custom barycenter verification test
add_executable(test_custom_barycenter test_custom_barycenter.cpp)
target_link_libraries(test_custom_barycenter PRIVATE delaunay_interfaces)

add_test(NAME CustomBarycenterTest COMMAND test_custom_barycenter)

# You can add more sophisticated testing with Google Test or Catch2
# For now, this is a minimal test setup
