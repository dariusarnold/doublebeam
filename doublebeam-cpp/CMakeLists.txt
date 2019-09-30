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
#find_package(xsimd)
#add_definitions(-DXTENSOR_ENABLE_XSIMD)

add_library(doublebeam_lib src/model.cpp src/utils.cpp src/raytracing.cpp src/twopoint.cpp
        src/raytracing_helpers.cpp src/raytracing_helpers.hpp src/ray.cpp src/ray.hpp
        src/beam.cpp src/beam.hpp)
target_link_libraries(doublebeam_lib PRIVATE Boost::boost stdc++fs)
target_link_libraries(doublebeam_lib PUBLIC xtensor ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_include_directories(doublebeam_lib PUBLIC src)
target_compile_options(doublebeam_lib PRIVATE -Werror -Wall -Wextra)
#target_compile_options(doublebeam_lib PUBLIC -ffast-math -march=native)

add_executable(doublebeam main.cpp)
target_compile_options(doublebeam PRIVATE -Werror -Wall -Wextra)
#target_link_libraries(doublebeam PRIVATE timing)
target_include_directories(doublebeam PRIVATE ${xtensor_blas_INCLUDE_DIRS})

add_subdirectory(src/timing)
target_link_libraries(doublebeam PRIVATE doublebeam_lib Boost::boost)

add_subdirectory(tests)

add_executable(scratch_1 /home/darius/.CLion2019.2/config/scratches/scratch_1.cpp)
target_link_libraries(scratch_1 PRIVATE timing)
add_executable(scratch_2 /home/darius/.CLion2019.2/config/scratches/scratch_2.cpp)
target_link_libraries(scratch_2 xtensor)
