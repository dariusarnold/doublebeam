project(${CMAKE_PROJECT_NAME}_tests)

find_package(GTest)

add_library(UtilitiesForTesting testing_utils.hpp testing_utils.cpp)
# add current folder as include directory so tests in folder unit_tests dont have to use
# .. to include header.
target_include_directories(UtilitiesForTesting PUBLIC .)
target_link_libraries(UtilitiesForTesting PRIVATE GTest::GTest doublebeam_lib)

add_subdirectory(unit_tests)