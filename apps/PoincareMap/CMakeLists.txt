add_executable(PoincareMap main.cpp)

target_link_libraries(PoincareMap PRIVATE
    my_lib
    Boost::boost
)

target_compile_definitions(PoincareMap PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS PoincareMap DESTINATION bin)
