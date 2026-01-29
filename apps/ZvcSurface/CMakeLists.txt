add_executable(ZvcSurface main.cpp)

target_link_libraries(ZvcSurface PRIVATE my_lib)

target_compile_definitions(ZvcSurface PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(ZvcSurface PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

install(TARGETS ZvcSurface DESTINATION bin)
