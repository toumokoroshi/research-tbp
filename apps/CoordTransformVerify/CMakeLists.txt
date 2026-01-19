add_executable(CoordTransformVerify main.cpp)

target_link_libraries(CoordTransformVerify PRIVATE my_lib)

target_compile_definitions(CoordTransformVerify PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(CoordTransformVerify PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

install(TARGETS CoordTransformVerify DESTINATION bin)
