cmake_minimum_required(VERSION 3.5)
project(TEST VERSION 1.0 LANGUAGES CXX)

#=============================================================== Setup CMake
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_FLAGS --coverage)
# set(CMAKE_EXE_LINKER_FLAGS --coverage)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include_directories(
  ${CMAKE_SOURCE_DIR}/../algorithm
)

enable_testing()

add_executable(
  ${PROJECT_NAME}
  ../algorithm/DNA_Analyzer.cc
  ../algorithm/DNA_Analyzer.h
  tests.cc
)

configure_file(${CMAKE_SOURCE_DIR}/HIV-1_AF033819.3.txt ${CMAKE_BINARY_DIR}/HIV-1_AF033819.3.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/test_1.txt ${CMAKE_BINARY_DIR}/test_1.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/test_2.txt ${CMAKE_BINARY_DIR}/test_2.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/test_3.txt ${CMAKE_BINARY_DIR}/test_3.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/test_4.txt ${CMAKE_BINARY_DIR}/test_4.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/test_5.txt ${CMAKE_BINARY_DIR}/test_5.txt COPYONLY)

target_link_libraries(${PROJECT_NAME}
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})