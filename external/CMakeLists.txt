project(${CMAKE_PROJECT_NAME}_external)

# clone cnpy
set(CNPY_DIR ${CMAKE_SOURCE_DIR}/external/cnpy)
if (NOT EXISTS ${CNPY_DIR})
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} clone https://github.com/rogersce/cnpy.git ${CNPY_DIR})
endif()

add_subdirectory(${CNPY_DIR})
target_include_directories(cnpy SYSTEM INTERFACE cnpy)
