cmake_minimum_required(VERSION 3.10)
project(doublebeam_cpp CXX)

set(CMAKE_CXX_STANDARD 17)

# set boost path with BOOST_ROOT if not globally installed
find_package(Boost REQUIRED)
# if you have installed these libraries using anaconda set CMAKE_PREFIX_PATH to
# /path/to/your/condaenv/lib/cmake/
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)

add_library(doublebeam_lib src/model.cpp src/utils.cpp src/raytracing.cpp src/twopoint.cpp)
target_link_libraries(doublebeam_lib PRIVATE Boost::boost stdc++fs)
target_link_libraries(doublebeam_lib PUBLIC xtensor)
target_include_directories(doublebeam_lib PUBLIC src)
target_compile_options(doublebeam_lib PRIVATE -Werror -Wall -Wextra)

add_executable(doublebeam main.cpp)
target_compile_options(doublebeam PRIVATE -Werror -Wall -Wextra)
target_link_libraries(doublebeam PRIVATE timing)

add_subdirectory(src/timing)
target_link_libraries(doublebeam PRIVATE doublebeam_lib Boost::boost)

add_subdirectory(tests)
