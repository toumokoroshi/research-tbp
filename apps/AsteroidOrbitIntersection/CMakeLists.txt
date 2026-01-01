add_executable(AsteroidOrbitIntersection main.cpp)

target_link_libraries(AsteroidOrbitIntersection PRIVATE
    my_lib
)

# OpenMPの設定
if(OpenMP_CXX_FOUND)
    # コンパイル時にのみOpenMPフラグを追加
    target_compile_options(AsteroidOrbitIntersection PRIVATE ${OpenMP_CXX_FLAGS})

    # リンク時はライブラリのみを追加（フラグは追加しない）
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # Intel OpenMPライブラリを直接リンク
        target_link_libraries(AsteroidOrbitIntersection PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        # 他の環境では通常のOpenMPリンク
        target_link_libraries(AsteroidOrbitIntersection PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(AsteroidOrbitIntersection PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS AsteroidOrbitIntersection DESTINATION bin)
