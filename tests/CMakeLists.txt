project(${CMAKE_PROJECT_NAME}_tests)

# clone google test
set(GOOGLETEST_DIR ${CMAKE_SOURCE_DIR}/tests/lib/googletest)
if (NOT EXISTS ${GOOGLETEST_DIR})
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} clone https://github.com/google/googletest ${GOOGLETEST_DIR})
endif()

add_library(UtilitiesForTesting testing_utils.hpp testing_utils.cpp)
# add current folder as include directory so tests in folder unit_tests dont have to use
# .. to include header.
target_include_directories(UtilitiesForTesting PUBLIC .)
target_link_libraries(UtilitiesForTesting PRIVATE gtest doublebeam_lib)

add_subdirectory(lib/googletest)
add_subdirectory(unit_tests)