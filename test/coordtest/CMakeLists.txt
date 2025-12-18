add_executable(coordtest main.cpp)
add_executable(test_v2 test_v2.cpp)
add_executable(test_v2_realdata test_v2_realdata.cpp)

target_link_libraries(coordtest PRIVATE my_lib)
target_link_libraries(test_v2 PRIVATE my_lib)
target_link_libraries(test_v2_realdata PRIVATE my_lib)

# OpenMPの設定
if(OpenMP_CXX_FOUND)
    # コンパイル時にのみOpenMPフラグを追加
    target_compile_options(coordtest PRIVATE ${OpenMP_CXX_FLAGS})

    # リンク時はライブラリのみを追加（フラグは追加しない）
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # Intel OpenMPライブラリを直接リンク
        target_link_libraries(coordtest PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        # 他の環境では通常のOpenMPリンク
        target_link_libraries(coordtest PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(coordtest PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/test/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/test/data"
    PROJECT_BASE_DIR="${CMAKE_SOURCE_DIR}"
)

# set_target_properties(coordtest PROPERTIES OUTPUT_NAME "coordtest_v3")

install(TARGETS coordtest DESTINATION bin)
