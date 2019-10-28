#include <filesystem>

#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include "seismodata.hpp"
#include "printing.hpp"

class TestProjectLoading : public testing::TestWithParam<std::string> {
protected:
    std::filesystem::path project_dir = current_source_path(__FILE__) / "data" / GetParam();
};

TEST_P(TestProjectLoading, TestIfLoadsWithoutThrow) {
    ASSERT_NO_THROW(SeismoData s(project_dir));
}

TEST_P(TestProjectLoading, TestCorrectnessOfLoadedProject) {
    SeismoData s(project_dir);
    ASSERT_EQ(s.receivers().size(), 3);
    ASSERT_EQ(s.sources().size(), 2);
    // loading source and receiver file is already tested.
}

TEST_P(TestProjectLoading, TestRetrievingSeismograms) {
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

INSTANTIATE_TEST_SUITE_P(TestNormalTextFileLoading, TestProjectLoading,
                         testing::Values("sample_project_dir"));

INSTANTIATE_TEST_SUITE_P(TestBinaryFileLoading, TestProjectLoading,
                         testing::Values("sample_project_dir_binary"));


struct TestSeismogramCutData {
    Seismogram expected;
    double t0, t1;
};

std::ostream& operator<<(std::ostream& os, const Seismogram& seismogram) {
    os << "Seismogram(t = " << seismogram.timesteps << ", x = " << seismogram.data << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const TestSeismogramCutData& data) {
    os << "TestSeismogramCutData(expected = " << data.expected << ", t0 = " << data.t0 << ", t1 = " << data.t1
       << ")";
    return os;
}

class TestSeismogramCut : public ::testing::TestWithParam<TestSeismogramCutData> {
protected:
    Seismogram in{{0, 1, 2, 3, 4, 5, 6}, {1, 2.5, 3, 5, 7, 8.1, 9}};
};

TEST_P(TestSeismogramCut, CompareResultWithExpected) {
    auto [expected, t0, t1] = GetParam();
    auto out = cut(in, t0, t1);
    ASSERT_EQ(out.data.size(), expected.data.size()) << "Wrong size for cut seismogram.";
    EXPECT_EQ(out.data, expected.data);
}

INSTANTIATE_TEST_SUITE_P(TestEmptySeismogram, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{{}, 10, 11}));

INSTANTIATE_TEST_SUITE_P(TestKeepingFullSeismogram, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{
                             Seismogram{{0, 1, 2, 3, 4, 5, 6}, {1, 2.5, 3, 5, 7, 8.1, 9}}, -0.5,
                             10}));

INSTANTIATE_TEST_SUITE_P(TestKeepingFirstValue, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{Seismogram{{0}, {1}}, -0.5, 0.5}));

INSTANTIATE_TEST_SUITE_P(TestKeepingLastValue, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{Seismogram{{6}, {9}}, 5.9, 10}));

INSTANTIATE_TEST_SUITE_P(TestInclusivityOfT0, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{Seismogram{{1, 2}, {2.5, 3}}, 1,
                                                               2.5}));

INSTANTIATE_TEST_SUITE_P(TestKeepingInclusivityT1, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{Seismogram{{3, 4, 5}, {5, 7, 8.1}},
                                                               2.5, 5}));