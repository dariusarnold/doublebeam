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
    for (const auto& source : s.sources()) {
        for (const auto& receiver : s.receivers()) {
            ASSERT_EQ(s(source, receiver).size(), 1000);
        }
    }
}
