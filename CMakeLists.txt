cmake_minimum_required(VERSION 3.16)

project(BooleanFunctions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(BF_LIB_HEADERS
    src/BF.h
    src/Utils.h
)

set(BF_LIB_SOURCES
    src/BF.cpp
    src/BFTransforms.cpp
    src/BFMetrics.cpp
    src/BFGeneration.cpp
    src/BFImprove.cpp
    src/Utils.cpp
)

add_library(bf_lib STATIC
    ${BF_LIB_HEADERS}
    ${BF_LIB_SOURCES}
)
source_group("Header Files" FILES ${BF_LIB_HEADERS})
source_group("Source Files" FILES ${BF_LIB_SOURCES})

target_include_directories(bf_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(bf_main src/main.cpp)
target_link_libraries(bf_main PRIVATE bf_lib)



include(CTest)
if(BUILD_TESTING)
    include(FetchContent)

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    fetchcontent_declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    fetchcontent_makeavailable(googletest)

    add_executable(bf_tests test/bf_tests.cpp)
    target_link_libraries(bf_tests PRIVATE bf_lib GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(bf_tests)
endif()
