add_executable(test main.cpp)

target_link_libraries(test PRIVATE
    Boost::boost
)


# OpenMPの設定
if(OpenMP_CXX_FOUND)
    # コンパイル時にのみOpenMPフラグを追加
    target_compile_options(test PRIVATE ${OpenMP_CXX_FLAGS})

    # リンク時はライブラリのみを追加（フラグは追加しない）
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # Intel OpenMPライブラリを直接リンク
        target_link_libraries(test PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        # 他の環境では通常のOpenMPリンク
        target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(test PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_BINARY_DIR}/results"
)

install(TARGETS test DESTINATION bin)