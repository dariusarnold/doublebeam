include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(Unit_Tests_run
        test_utils.cpp
        test_integration.cpp test_velocity_model.cpp test_ray_tracing.cpp test_twopoint.cpp test_dynamic_raytracing.cpp)

target_link_libraries(Unit_Tests_run gtest gtest_main)
target_link_libraries(Unit_Tests_run UtilitiesForTesting)
target_link_libraries(Unit_Tests_run doublebeam_lib Boost::boost)