add_executable(hello main.cpp)

target_link_libraries(hello PRIVATE
  my_lib
)

target_compile_definitions(hello PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS hello DESTINATION bin)