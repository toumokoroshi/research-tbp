cmake_minimum_required(VERSION 3.20)

# PeriodicOrbitContinuation アプリケーション
add_executable(PeriodicOrbitContinuation main.cpp)

# CONFIG_DIR と OUTPUT_DIR の定義
target_compile_definitions(PeriodicOrbitContinuation PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

# インクルードパス
target_include_directories(PeriodicOrbitContinuation PRIVATE ${CMAKE_SOURCE_DIR}/lib/include)

# toml++ の取得と設定
include(FetchContent)
FetchContent_Declare(
  tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
  GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

# リンク
target_link_libraries(PeriodicOrbitContinuation PRIVATE tomlplusplus::tomlplusplus)

# C++20を使用（std::numbersに必要）
target_compile_features(PeriodicOrbitContinuation PRIVATE cxx_std_20)

# Boost.Odeintが必要な場合（rtbp.hppが使用）
find_package(Boost REQUIRED)
target_link_libraries(PeriodicOrbitContinuation PRIVATE Boost::boost)

# 出力先を bin/ に設定
set_target_properties(PeriodicOrbitContinuation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(PeriodicOrbitContinuation PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled for PeriodicOrbitContinuation")
endif()

# Windowsの場合、マルチプロセッサコンパイルを有効化
if(MSVC)
    target_compile_options(PeriodicOrbitContinuation PRIVATE /MP)
endif()
