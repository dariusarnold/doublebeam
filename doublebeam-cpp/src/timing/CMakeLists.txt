cmake_minimum_required(VERSION 3.10)
project(timing)
set(CMAKE_CXX_STANDARD 17)

add_library(timing timing.cpp)
target_compile_options(timing PRIVATE -Werror -Wall -Wextra)
target_include_directories(timing PUBLIC .)