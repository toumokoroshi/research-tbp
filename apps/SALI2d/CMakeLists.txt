add_executable(SALI2d main.cpp)
target_link_libraries(SALI2d PRIVATE
    mylib

    # Eigen3::Eigen
    # Boost::boost
    # OpenMP::OpenMP_CXX
)
target_compile_definitions(SALI2d PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/output"
)

install(TARGETS SALI2d DESTINATION bin)