# Benchmarks (not tests: they print timings, assert nothing; run manually
# or via scripts/run_benchmarks.sh).

add_executable(benchmark_speed benchmark_speed.cpp)
target_link_libraries(benchmark_speed PRIVATE delaunay_interfaces)
