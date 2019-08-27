cmake_minimum_required(VERSION 3.14)
project(doublebeam_cpp)

set(CMAKE_CXX_STANDARD 17)

add_library(timing src/timing/timing.cpp)

add_executable(doublebeam main.cpp)

target_link_libraries(doublebeam PRIVATE timing)