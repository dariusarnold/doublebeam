#include <filesystem>

#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include "seismodata.hpp"

class TestProjectLoading : public testing::Test {
protected:
    std::filesystem::path project_dir =
        current_source_path(__FILE__) / "data" / "sample_project_dir";
};

TEST_F(TestProjectLoading, TestIfLoadsWithoutThrow) {
    ASSERT_NO_THROW(SeismoData s(project_dir));
}

TEST_F(TestProjectLoading, TestCorrectnessOfLoadedProject) {
    SeismoData s(project_dir);
    ASSERT_EQ(s.receivers().size(), 3);
    ASSERT_EQ(s.sources().size(), 2);
    // loading source and receiver file is already tested.
}

TEST_F(TestProjectLoading, TestRetrievingSeismograms) {
    SeismoData s(project_dir);
    // I replaced the first value of every seismogram with an increasing index to test if they are
    // loaded in the correct order.
    auto i = 1;
    for (const auto& source : s.sources()) {
        for (const auto& receiver : s.receivers()) {
            auto seismogram = s(source, receiver);
            ASSERT_EQ(seismogram.data.size(), 1000);
            EXPECT_EQ(seismogram.data[0], i) << "Loaded wrong seismogram for source " << source
                                             << ", receiver " << receiver << ".";
            ++i;
        }
    }
}
