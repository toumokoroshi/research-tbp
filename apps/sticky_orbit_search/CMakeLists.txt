
set(TARGET_NAME sticky_orbit_search)

add_executable(${TARGET_NAME}
    main.cpp
)

target_link_libraries(${TARGET_NAME} PRIVATE
    my_lib
)

# OpenMPの設定
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_compile_options(${TARGET_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
    if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        target_link_libraries(${TARGET_NAME} PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        target_link_libraries(${TARGET_NAME} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(${TARGET_NAME} PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(${TARGET_NAME} PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)
