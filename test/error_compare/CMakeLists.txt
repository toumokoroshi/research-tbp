add_executable(error_compare main.cpp)

target_link_libraries(error_compare PRIVATE
  my_lib
)

set_target_properties(error_compare PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

target_compile_definitions(error_compare PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS error_compare DESTINATION bin)
