add_executable(SALI3d main.cpp)
target_link_libraries(SALI3d PRIVATE
    mylib

    # Eigen3::Eigen
    Boost::boost
    OpenMP::OpenMP_CXX
)
target_compile_definitions(SALI3d PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_BINARY_DIR}/results"
)

install(TARGETS SALI3d DESTINATION bin)