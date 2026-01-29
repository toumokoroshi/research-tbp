cmake_minimum_required(VERSION 3.20)

# StableManifoldComputation アプリケーション
add_executable(StableManifoldComputation main.cpp)

# CONFIG_DIR と OUTPUT_DIR の定義
target_compile_definitions(StableManifoldComputation PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
    BOOST_NO_EXCEPTIONS
    BOOST_EXCEPTION_DISABLE
)

# インクルードパス
target_include_directories(StableManifoldComputation PRIVATE 
    ${CMAKE_SOURCE_DIR}/lib/include
    ${eigen_SOURCE_DIR}
)

# リンク
target_link_libraries(StableManifoldComputation PRIVATE 
    tomlplusplus::tomlplusplus
    Boost::boost
)

# C++20を使用
target_compile_features(StableManifoldComputation PRIVATE cxx_std_20)

# 出力先を bin/ に設定
set_target_properties(StableManifoldComputation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(StableManifoldComputation PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled for StableManifoldComputation")
endif()

# Windowsの場合、マルチプロセッサコンパイルを有効化
if(MSVC)
    target_compile_options(StableManifoldComputation PRIVATE /MP)
endif()
