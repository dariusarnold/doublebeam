cmake_minimum_required(VERSION 3.14)
project(doublebeam_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/findFFTW")

# set boost path with BOOST_ROOT if not globally installed
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt 5 REQUIRED)
find_package(FFTW COMPONENTS FLOAT_LIB DOUBLE_LIB REQUIRED)

message("-- Setting up nanoflann.")
execute_process(COMMAND git submodule update --init -- external/nanoflann
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(external/nanoflann EXCLUDE_FROM_ALL)

set(WARNINGS_GCC -Wall -Wextra -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wduplicated-cond -Wuseless-cast -Wdouble-promotion -Wduplicated-branches -Wmisleading-indentation
        -Wshadow -Wextra-semi)

set(WARNINGS_CLANG -Wall -Wextra -Wnon-virtual-dtor -Wold-style-cast -Wunused -Woverloaded-virtual
        -Wdouble-promotion -Woverloaded-shift-op-parentheses -Wshadow -Wextra-semi -Wno-c++98-compat-extra-semi)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(WARNINGS ${WARNINGS_CLANG})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(WARNINGS ${WARNINGS_GCC})
endif()

add_library(doublebeam_lib src/model.cpp src/utils.cpp src/raytracing.cpp src/twopoint.cpp
        src/raytracing_helpers.cpp src/ray.cpp src/beam.cpp src/printing.cpp src/io.cpp
        src/seismodata.cpp src/fft.cpp src/eigen_helpers.cpp src/strong_types.hpp src/units.hpp
        src/raytracing_types.cpp src/kdtree.cpp src/kdtree.hpp)
target_link_libraries(doublebeam_lib PRIVATE stdc++fs FFTW::Float FFTW::Double fmt::fmt)
target_link_libraries(doublebeam_lib PUBLIC Eigen3::Eigen cnpy nanoflann::nanoflann Boost::boost)
target_include_directories(doublebeam_lib PUBLIC src)
target_compile_options(doublebeam_lib PRIVATE ${WARNINGS})
#target_compile_options(doublebeam_lib PUBLIC -ffast-math -march=native)

add_executable(doublebeam main.cpp src/doublebeam.cpp src/doublebeam.hpp)
target_compile_options(doublebeam PRIVATE ${WARNINGS})
target_link_libraries(doublebeam PRIVATE doublebeam_lib Boost::boost Boost::program_options Eigen3::Eigen cnpy)

add_executable(convert_to_binary convert_to_binary.cpp)
target_link_libraries(convert_to_binary stdc++fs doublebeam_lib)

# Set ENABLE_LTO to false to disable link time optimization, disabling the compiler
# from optimizing across translation units. This results in lower linking times,
# but might increase runtime of the program.
option(ENABLE_LTO "Set to false to disable link time optimization." true)
if (ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error LANGUAGES CXX)
endif()
if (ENABLE_LTO AND supported)
    message(STATUS "IPO/LTO enabled")
    set_property(TARGET doublebeam_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET doublebeam PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_link_options(doublebeam_lib PUBLIC -flto)
    endif()
else ()
    message(STATUS "IPO/LTO not supported; <${error}>")
endif ()


# Basic timing library runs one function and gives average runtime.
#add_subdirectory(src/timing)
# load vendored in external libraries not installed globally by the package manager
add_subdirectory(external)
# load unit tests
add_subdirectory(tests)


# Check if all headers are standalone includable.
# Currently this is only tested when cmake is called, since I haven't found a way to export this as a target.
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_QUIET ON)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_REQUIRED_DEFINITIONS -DBOOST_ROOT=/opt/boost_1_71/ -DCMAKE_PREFIX_PATH=/home/darius/.anaconda3/lib/cmake/)
set(CMAKE_REQUIRED_LIBRARIES stdc++fs Boost::boost Eigen3::Eigen FFTW::Float FFTW::Double cnpy)
file(GLOB HEADERS src/*.hpp)
if (FALSE)
    foreach (HEADER ${HEADERS})
        get_filename_component(HEADERNAME ${HEADER} NAME)
        message(STATUS "Standalone compilation of " ${HEADERNAME})
        check_include_file_cxx(${HEADER} DID_HEADER_BUILD)
        if (NOT DID_HEADER_BUILD)
            message(WARNING "Failed to build " ${HEADERNAME} " as a standalone header.\n"
                    "Check '" ${CMAKE_BINARY_DIR} "/CMakeFiles/CMakeError.log' for a line starting with\n"
                    "Determining if the include file " ${HEADER})
        endif ()
        # The result of check_include_file_cxx is saved in a cache variable. If the variable is already defined, for example
        # after the first loop iteration, the whole check is skipped. Because of that we reset the variable.
        # see https://gitlab.kitware.com/cmake/cmake/issues/18492
        # and line 46 in https://gitlab.kitware.com/cmake/cmake/blob/v3.14.6/Modules/CheckIncludeFileCXX.cmake
        unset(DID_HEADER_BUILD CACHE)
    endforeach ()
endif ()
