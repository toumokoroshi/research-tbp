add_executable(ChaosIndicator main.cpp)

target_link_libraries(ChaosIndicator PRIVATE my_lib)

# OpenMP setup (same pattern as TrajectorySali)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_compile_options(ChaosIndicator PRIVATE ${OpenMP_CXX_FLAGS})
    
    # Intel OpenMP library handling for MSVC
    if(MSVC AND OpenMP_iomp5md_LIBRARY)
        target_link_libraries(ChaosIndicator PRIVATE ${OpenMP_iomp5md_LIBRARY})
    else()
        target_link_libraries(ChaosIndicator PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_definitions(ChaosIndicator PRIVATE
    CONFIG_DIR="${CMAKE_SOURCE_DIR}/configs"
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/data"
)

install(TARGETS ChaosIndicator DESTINATION bin)
