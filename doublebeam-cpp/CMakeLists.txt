cmake_minimum_required(VERSION 3.10)
project(doublebeam_cpp CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

add_executable(doublebeam main.cpp src/model.cpp src/utils.cpp)
target_compile_options(doublebeam PRIVATE -Werror -Wall -Wextra -Wpedantic)
add_subdirectory(src/timing)
target_link_libraries(doublebeam PRIVATE Boost::boost timing)