add_executable(Unit_Tests_run
        test_utils.cpp test_raytracing_helpers.cpp test_velocity_model.cpp test_ray_tracing.cpp test_twopoint.cpp
        test_dynamic_raytracing.cpp test_printing.cpp test_io.cpp test_seismodata.cpp test_fft.cpp data/gtest.cpp
        test_eigen_helpers.cpp)

target_link_libraries(Unit_Tests_run GTest::GTest UtilitiesForTesting doublebeam_lib Boost::boost)

# copy data for tests into the build folder so relative paths can be used
add_custom_command(TARGET Unit_Tests_run POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/
        $<TARGET_FILE_DIR:Unit_Tests_run>/data)