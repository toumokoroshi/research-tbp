cmake_minimum_required(VERSION 3.15)

project(TestHaloBifurcation CXX)

add_executable(test_halo_bifurcation test_halo_bifurcation.cpp)

target_include_directories(test_halo_bifurcation PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/include
)

target_link_libraries(test_halo_bifurcation PRIVATE
  Boost::boost
  tomlplusplus::tomlplusplus
  Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(test_halo_bifurcation PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME test_halo_bifurcation COMMAND test_halo_bifurcation)

# L2 Lyapunov継続と分岐検出テスト
add_executable(test_l2_continuation test_l2_continuation.cpp)

target_include_directories(test_l2_continuation PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/include
)

target_link_libraries(test_l2_continuation PRIVATE
  Boost::boost
  tomlplusplus::tomlplusplus
  Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(test_l2_continuation PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME test_l2_continuation COMMAND test_l2_continuation)

# L2 既知初期条件検証テスト
add_executable(test_l2_known_ic test_l2_known_ic.cpp)

target_include_directories(test_l2_known_ic PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/include
)

target_link_libraries(test_l2_known_ic PRIVATE
  Boost::boost
  tomlplusplus::tomlplusplus
  Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(test_l2_known_ic PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME test_l2_known_ic COMMAND test_l2_known_ic)

# Halo軌道継続テスト（分岐点から生成）
add_executable(test_halo_continuation test_halo_continuation.cpp)

target_include_directories(test_halo_continuation PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/include
)

target_link_libraries(test_halo_continuation PRIVATE
  Boost::boost
  tomlplusplus::tomlplusplus
  Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(test_halo_continuation PRIVATE OpenMP::OpenMP_CXX)
endif()

add_test(NAME test_halo_continuation COMMAND test_halo_continuation)

