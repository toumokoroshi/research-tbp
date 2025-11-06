add_executable(SSBI2SEBR main.cpp)

target_link_libraries(SSBI2SEBR PRIVATE
  my_lib
)

target_compile_definitions(SSBI2SEBR PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS SSBI2SEBR DESTINATION bin)
