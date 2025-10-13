cmake_minimum_required(VERSION 3.15)

# Windows環境でIntel oneAPIを使用する場合の設定
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # MSVCプラットフォームファイルを使わないように設定
    set(CMAKE_CXX_COMPILER_FORCED TRUE)
    set(CMAKE_C_COMPILER_FORCED TRUE)
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

# # Intel MKLの設定
# set(MKL_ROOT "${ONEAPI_ROOT}/mkl/latest")
# set(MKL_INTERFACE "lp64")  # または "ilp64" (64bit整数)
# set(MKL_THREADING "tbb")   # または "sequential", "intel_thread"
# set(MKL_MPI "openmpi")     # MPIを使う場合

# # Intel MKLの検索
# find_package(MKL CONFIG PATHS ${MKL_ROOT})

# if(MKL_FOUND)
#     message(STATUS "Intel MKL found")
#     set(BLAS_LIBRARIES MKL::MKL)
# else()
#     message(STATUS "Intel MKL not found, using system BLAS")
#     find_package(BLAS REQUIRED)
# endif()

# # Intel TBBの検索
# find_package(TBB CONFIG PATHS "${ONEAPI_ROOT}/tbb/latest")

# if(TBB_FOUND)
#     message(STATUS "Intel TBB found")
# endif()

# コンパイラフラグ
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -qopt-report=5 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wshadow -Wconversion -Wtrampolines -Wimplicit-fallthrough -Wno-dev")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# -------------------------------------------------------------
# 外部ライブラリの検索
# -------------------------------------------------------------
# # Eigen (ヘッダーオンリーなので見つけるだけでOK)
# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Boost
find_package(Boost 1.71 REQUIRED)
# Intel oneAPI OpenMP設定
if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(ONEAPI_ROOT "C:/Program Files (x86)/Intel/oneAPI")
    
    # compiler/latestディレクトリを使用
    if(EXISTS "${ONEAPI_ROOT}/compiler/latest")
        set(INTEL_COMPILER_DIR "${ONEAPI_ROOT}/compiler/latest")
        
        # OpenMPフラグとライブラリを設定
        if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
            # 新しいIntel LLVMベースのコンパイラ (icx/icpx)
            set(OpenMP_CXX_FLAGS "-fiopenmp")
            set(OpenMP_CXX_LIB_NAMES "iomp5md")
            set(OpenMP_iomp5md_LIBRARY "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
        else()
            message(WARNING "Unexpected compiler ID: ${CMAKE_CXX_COMPILER_ID}")
        endif()
        
        # ライブラリの存在確認
        if(EXISTS "${OpenMP_iomp5md_LIBRARY}")
            message(STATUS "Found Intel OpenMP library: ${OpenMP_iomp5md_LIBRARY}")
        else()
            message(WARNING "Intel OpenMP library not found at expected location")
            # 32ビット版も試す
            set(OpenMP_iomp5md_LIBRARY "${INTEL_COMPILER_DIR}/compiler/lib/ia32_win/libiomp5md.lib")
        endif()
    endif()
endif()

find_package(OpenMP REQUIRED)
 
# -------------------------------------------------------------
# サブディレクトリの追加
# -------------------------------------------------------------
add_subdirectory(lib)
add_subdirectory(apps)