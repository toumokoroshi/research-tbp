add_executable(TrajectorySali main.cpp)

target_link_libraries(TrajectorySali PRIVATE my_lib)

# OpenMPの設定
if(OpenMP_CXX_FOUND)
    # コンパイル時にのみOpenMPフラグを追加
    target_compile_options(TrajectorySali PRIVATE ${OpenMP_CXX_FLAGS})

    # リンク時はライブラリのみを追加（フラグは追加しない）
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # Intel OpenMPライブラリを直接リンク
        target_link_libraries(TrajectorySali PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        # 他の環境では通常のOpenMPリンク
        target_link_libraries(TrajectorySali PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(TrajectorySali PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

# set_target_properties(TrajectorySali PROPERTIES OUTPUT_NAME "TrajectorySali_v3")

install(TARGETS TrajectorySali DESTINATION bin)
