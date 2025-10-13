cmake_minimum_required(VERSION 3.15)

# WindowsでIntel LLVMコンパイラを使う場合、MSVCプラットフォームとして扱う
if(WIN32 AND (CMAKE_CXX_COMPILER MATCHES "icpx" OR CMAKE_CXX_COMPILER MATCHES "icx"))
    set(CMAKE_CXX_COMPILER_ID "MSVC" CACHE STRING "" FORCE)
    set(CMAKE_C_COMPILER_ID "MSVC" CACHE STRING "" FORCE)
endif()

# プロジェクト名
project(MySimulation VERSION 1.0.0 LANGUAGES CXX)

# C++の標準を設定 (C++17を推奨)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Intel oneAPIの検出
if(DEFINED ENV{ONEAPI_ROOT})
    set(ONEAPI_ROOT $ENV{ONEAPI_ROOT})
    message(STATUS "Intel oneAPI detected at: ${ONEAPI_ROOT}")
else()
    set(ONEAPI_ROOT "/opt/intel/oneapi")  # デフォルトパス
endif()

# WindowsでIntel LLVMコンパイラを使用する場合の設定
if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Configuring for Intel LLVM compiler on Windows (MSVC-compatible mode)")
    
    # MSVC互換モードのフラグを使用
    set(CMAKE_CXX_FLAGS_RELEASE "/O3 /Qipo /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /Wall")
    
    # ランタイムライブラリの設定
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
    
else()
    # Unix系またはGNU/Clangモードの場合
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wshadow -Wconversion")
endif()

# ビルドタイプが指定されていない場合はDebugに設定
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# -------------------------------------------------------------
# 外部ライブラリの検索
# -------------------------------------------------------------
# Boost
find_package(Boost 1.71 REQUIRED)

# Intel oneAPI OpenMP設定
if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    set(ONEAPI_ROOT "C:/Program Files (x86)/Intel/oneAPI")
    
    if(EXISTS "${ONEAPI_ROOT}/compiler/latest")
        set(INTEL_COMPILER_DIR "${ONEAPI_ROOT}/compiler/latest")
        
        # OpenMPフラグとライブラリを設定 (MSVC互換形式)
        set(OpenMP_CXX_FLAGS "/Qopenmp")
        set(OpenMP_CXX_LIB_NAMES "iomp5md")
        set(OpenMP_iomp5md_LIBRARY "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
        
        # ライブラリの存在確認
        if(EXISTS "${OpenMP_iomp5md_LIBRARY}")
            message(STATUS "Found Intel OpenMP library: ${OpenMP_iomp5md_LIBRARY}")
        else()
            message(WARNING "Intel OpenMP library not found at expected location")
        endif()
    endif()
endif()

find_package(OpenMP REQUIRED)
 
# -------------------------------------------------------------
# サブディレクトリの追加
# -------------------------------------------------------------
add_subdirectory(lib)
add_subdirectory(apps)