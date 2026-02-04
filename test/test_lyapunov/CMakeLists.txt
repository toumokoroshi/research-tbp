# テストアプリケーション

# Lyapunov初期条件テスト
add_executable(test_lyapunov test_lyapunov.cpp)

target_compile_features(test_lyapunov PRIVATE cxx_std_20)

# my_lib を使用（Boost, Eigen, インクルードパスを含む）
target_link_libraries(test_lyapunov PRIVATE my_lib)

set_target_properties(test_lyapunov PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

if(MSVC)
    target_compile_options(test_lyapunov PRIVATE /MP)
endif()
