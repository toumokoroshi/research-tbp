add_executable(JacobiIntegralCalc main.cpp)

target_link_libraries(JacobiIntegralCalc PRIVATE my_lib)

target_compile_definitions(JacobiIntegralCalc PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(JacobiIntegralCalc PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

install(TARGETS JacobiIntegralCalc DESTINATION bin)
