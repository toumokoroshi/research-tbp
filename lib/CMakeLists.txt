# 共有ライブラリのソースファイル
set(LIB_SOURCES
    src/crtbp.cpp
    src/rtbp.cpp
    # src/solver.cpp
    # src/io_utils.cpp
)

# 共有ライブラリの作成
add_library(mylib STATIC ${LIB_SOURCES})

# インクルードディレクトリの設定
target_include_directories(mylib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# C++標準の設定(一応全体で統一してみる)
# target_compile_features(mylib PUBLIC cxx_std_17)
set_property(TARGET mylib PROPERTY CXX_STANDARD 17)
set_property(TARGET mylib PROPERTY CXX_STANDARD_REQUIRED ON)

# 外部ライブラリとのリンク
target_link_libraries(mylib PUBLIC Boost::headers)
target_link_libraries(mylib PUBLIC OpenMP::OpenMP_CXX)
# endif()
# if (Boost_FOUND)
#     target_include_directories(mylib PUBLIC ${Boost_INCLUDE_DIRS})
#     target_link_libraries(mylib PUBLIC ${Boost_LIBRARIES})

# # mylibをインストールする設定（オプション）
# # これにより、ビルドしたライブラリとヘッダーをシステムにインストールできる
# install(TARGETS mylib DESTINATION lib)
# install(DIRECTORY include/mylib DESTINATION include)