# C++の標準規格を指定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Boost REQUIRED)

# ヘッダオンリーライブラリを INTERFACE ターゲットとして定義する
add_library(my_lib INTERFACE)

# my_header_lib にインクルードパスを設定する
target_include_directories(my_lib
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        $<INSTALL_INTERFACE:include>
)

#my_header_lib にBoostへの依存関係を設定する
target_link_libraries(my_lib
    INTERFACE
        Boost::headers
        tomlplusplus::tomlplusplus
)
