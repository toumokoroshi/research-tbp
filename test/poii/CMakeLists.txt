add_executable(poii main.cpp)

target_link_libraries(poii PRIVATE my_lib)

# OpenMPの設定
if(OpenMP_CXX_FOUND)
    # コンパイル時にのみOpenMPフラグを追加
    target_compile_options(poii PRIVATE ${OpenMP_CXX_FLAGS})

    # リンク時はライブラリのみを追加（フラグは追加しない）
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # Intel OpenMPライブラリを直接リンク
        target_link_libraries(poii PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        # 他の環境では通常のOpenMPリンク
        target_link_libraries(poii PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(poii PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

# set_target_properties(poii PROPERTIES OUTPUT_NAME "poii_v3")

install(TARGETS poii DESTINATION bin)
