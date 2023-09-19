cmake_minimum_required(VERSION 3.20)
project(LinpackTest)

set(CMAKE_CXX_STANDARD 11)

add_executable(
        LinpackTest
        main.cpp
        Matrix.hpp
        Matrix.cpp
        LinpackTester.hpp
        SimpleIterationTester.hpp
        matrix_error.hpp
        SimpleIterationParallelTester.hpp
        Matrix.cpp
        LinpackResult.cpp
        LinpackResult.h)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
        MatrixTest
        MatrixTest.cpp
        Matrix.hpp
)

target_link_libraries(
        MatrixTest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(MatrixTest)




