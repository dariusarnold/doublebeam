cmake_minimum_required(VERSION 3.14)
project(doublebeam_cpp CXX)

set(CMAKE_CXX_STANDARD 17)

# set boost path with BOOST_ROOT if not globally installed
find_package(Boost REQUIRED)
# if you have installed these libraries using anaconda set CMAKE_PREFIX_PATH to
# /path/to/your/condaenv/lib/cmake/
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)
add_definitions(-DHAVE_CBLAS=1)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
# the following is taken from https://github.com/egpbos/findFFTW
configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

find_package(FFTW REQUIRED)
#find_package(xsimd)
#add_definitions(-DXTENSOR_ENABLE_XSIMD)

add_library(doublebeam_lib src/model.cpp src/utils.cpp src/raytracing.cpp src/twopoint.cpp
        src/raytracing_helpers.cpp src/raytracing_helpers.hpp src/ray.cpp src/ray.hpp
        src/beam.cpp src/beam.hpp src/printing.cpp src/io.cpp src/io.hpp src/seismodata.hpp src/seismodata.cpp)
target_link_libraries(doublebeam_lib PRIVATE Boost::boost stdc++fs)
target_link_libraries(doublebeam_lib PUBLIC xtensor ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_include_directories(doublebeam_lib PUBLIC src)
target_compile_options(doublebeam_lib PRIVATE -Werror -Wall -Wextra)
#target_compile_options(doublebeam_lib PUBLIC -ffast-math -march=native)

add_executable(doublebeam main.cpp src/doublebeam.cpp src/doublebeam.hpp)
target_compile_options(doublebeam PRIVATE -Werror -Wall -Wextra)
target_link_libraries(doublebeam PRIVATE fftw3)
target_include_directories(doublebeam PRIVATE ${xtensor_blas_INCLUDE_DIRS} timing)

add_subdirectory(src/timing)
target_link_libraries(doublebeam PRIVATE doublebeam_lib Boost::boost timing)
target_include_directories(doublebeam PRIVATE FFTW_INCLUDE_DIRS)

add_subdirectory(tests)