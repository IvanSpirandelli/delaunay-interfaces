add_executable(simple_example simple_example.cpp)
target_link_libraries(simple_example PRIVATE delaunay_interfaces)

# Optional: install examples
install(TARGETS simple_example
    RUNTIME DESTINATION bin/examples
)
