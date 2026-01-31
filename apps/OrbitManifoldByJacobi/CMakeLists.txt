add_executable(OrbitManifoldByJacobi main.cpp)

target_include_directories(OrbitManifoldByJacobi PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/include
)

target_compile_definitions(OrbitManifoldByJacobi PRIVATE
  CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
)

target_link_libraries(OrbitManifoldByJacobi PRIVATE
  Boost::boost
  tomlplusplus::tomlplusplus
  Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(OrbitManifoldByJacobi PRIVATE OpenMP::OpenMP_CXX)
endif()

set_target_properties(OrbitManifoldByJacobi PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
)
