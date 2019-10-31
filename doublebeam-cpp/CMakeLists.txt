cmake_minimum_required(VERSION 3.14)
project(doublebeam_cpp CXX)

set(CMAKE_CXX_STANDARD 17)

# set boost path with BOOST_ROOT if not globally installed
find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
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
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if (result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
else ()
    message("CMake step for findFFTW completed (${result}).")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download)
if (result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif ()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

find_package(FFTW REQUIRED)
#find_package(xsimd)
#add_definitions(-DXTENSOR_ENABLE_XSIMD)

set(WARNINGS -Werror -Wall -Wextra -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wduplicated-cond -Wuseless-cast -Wdouble-promotion -Wduplicated-branches -Wmisleading-indentation -Weffc++)

add_library(doublebeam_lib src/model.cpp src/utils.cpp src/raytracing.cpp src/twopoint.cpp
        src/raytracing_helpers.cpp src/raytracing_helpers.hpp src/ray.cpp src/ray.hpp
        src/beam.cpp src/beam.hpp src/printing.cpp src/io.cpp src/io.hpp src/seismodata.hpp src/seismodata.cpp src/fft.cpp src/fft.hpp)
target_link_libraries(doublebeam_lib PRIVATE Boost::boost stdc++fs fftw3)
target_link_libraries(doublebeam_lib PUBLIC xtensor ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_include_directories(doublebeam_lib PUBLIC src)
target_include_directories(doublebeam_lib PRIVATE FFTW_INCLUDE_DIRS)
target_compile_options(doublebeam_lib PRIVATE ${WARNINGS})
#target_compile_options(doublebeam_lib PUBLIC -ffast-math -march=native)

add_executable(doublebeam main.cpp src/doublebeam.cpp src/doublebeam.hpp)
target_compile_options(doublebeam PRIVATE ${WARNINGS})
target_include_directories(doublebeam PRIVATE ${xtensor_blas_INCLUDE_DIRS} timing)
target_link_libraries(doublebeam PRIVATE doublebeam_lib Boost::boost timing Eigen3::Eigen)
add_subdirectory(src/timing)

add_executable(convert_to_binary convert_to_binary.cpp)
target_link_libraries(convert_to_binary stdc++fs doublebeam_lib)

add_subdirectory(tests)


# Check if all headers are standalone includable.
# Currently this is only tested when cmake is called, since I haven't found a way to export this as a target.
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_QUIET ON)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_REQUIRED_DEFINITIONS -DBOOST_ROOT=/opt/boost_1_71/ -DCMAKE_PREFIX_PATH=/home/darius/.anaconda3/lib/cmake/)
set(CMAKE_REQUIRED_LIBRARIES stdc++fs Boost::boost Eigen3::Eigen fftw3 xtensor ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
file(GLOB HEADERS src/*.hpp)
message(${CMAKE_BINARY_DIR})
foreach(HEADER ${HEADERS})
    get_filename_component(HEADERNAME ${HEADER} NAME)
    message(STATUS "Standalone compilation of " ${HEADERNAME})
    check_include_file_cxx(${HEADER} DID_HEADER_BUILD)
    if (NOT DID_HEADER_BUILD)
            message(WARNING "Failed to build " ${HEADERNAME} " as a standalone header.\n"
                    "Check '" ${CMAKE_BINARY_DIR} "/CMakeFiles/CMakeError.log' for a line starting with\n"
                    "Determining if the include file " ${HEADER})
    endif()
    # The result of check_include_file_cxx is saved in a cache variable. If the variable is already defined, for example
    # after the first loop iteration, the whole check is skipped. Because of that we reset the variable.
    # see https://gitlab.kitware.com/cmake/cmake/issues/18492
    # and line 46 in https://gitlab.kitware.com/cmake/cmake/blob/v3.14.6/Modules/CheckIncludeFileCXX.cmake
    unset(DID_HEADER_BUILD CACHE)
endforeach()
