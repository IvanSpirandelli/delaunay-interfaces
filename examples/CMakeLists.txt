add_executable(cpp_api cpp_api.cpp)
target_link_libraries(cpp_api PRIVATE delaunay_interfaces)

# Optional: install examples
install(TARGETS cpp_api
    RUNTIME DESTINATION bin/examples
)
