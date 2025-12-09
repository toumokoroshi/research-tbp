add_executable(TrajectoryCalc main.cpp)

target_link_libraries(TrajectoryCalc PRIVATE
  my_lib
)

target_compile_definitions(TrajectoryCalc PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS TrajectoryCalc DESTINATION bin)
