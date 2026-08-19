# Benchmarks (not tests: they print timings, assert nothing; run manually
# or via scripts/run_benchmarks.sh).

add_executable(benchmark_speed benchmark_speed.cpp)
target_link_libraries(benchmark_speed PRIVATE delaunay_interfaces)

# Compares the production radix bucket sort (src/filtration_sort.hpp)
# against std::sort on real pipeline buckets; needs the internal headers.
add_executable(benchmark_radix_sort benchmark_radix_sort.cpp)
target_link_libraries(benchmark_radix_sort PRIVATE delaunay_interfaces)
target_include_directories(benchmark_radix_sort PRIVATE ${CMAKE_SOURCE_DIR}/src)
