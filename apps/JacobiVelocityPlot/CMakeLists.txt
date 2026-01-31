add_executable(JacobiVelocityPlot main.cpp)

target_link_libraries(JacobiVelocityPlot PRIVATE my_lib)

target_compile_definitions(JacobiVelocityPlot PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(JacobiVelocityPlot PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

install(TARGETS JacobiVelocityPlot DESTINATION bin)
