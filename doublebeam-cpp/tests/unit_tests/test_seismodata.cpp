#include <filesystem>

#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include "seismodata.hpp"

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
    Seismogram in;
    Seismogram expected;
    std::vector<double> t;
    double t0, t1;
};

class TestSeismogramCut : public ::testing::TestWithParam<TestSeismogramCutData> {
protected:
};

TEST_P(TestSeismogramCut, CompareResultWithExpected) {
    auto [in, expected, t, t0, t1] = GetParam();
    auto out = cut(in, t, t0, t1);
    ASSERT_EQ(out.data.size(), expected.data.size()) << "Wrong size for cut seismogram.";
    EXPECT_EQ(out.data, expected.data);
}

INSTANTIATE_TEST_SUITE_P(TestEmptySeismogram, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{{}, {}, {}, 0, 1}));

INSTANTIATE_TEST_SUITE_P(
    TestKeepingFullSeismogram, TestSeismogramCut,
    testing::Values(TestSeismogramCutData{
        {{0, 1, 2, 3, 4, 5, 6}}, {{0, 1, 2, 3, 4, 5, 6}}, {1, 2.5, 3, 5, 7, 8.1, 9}, 0.5, 10}));

INSTANTIATE_TEST_SUITE_P(TestKeepingFirstValue, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{
                             {{0, 1, 2, 3, 4, 5, 6}}, {{0}}, {1, 2.5, 3, 5, 7, 8.1, 9}, 0.5, 1}));

INSTANTIATE_TEST_SUITE_P(TestKeepingLastValue, TestSeismogramCut,
                         testing::Values(TestSeismogramCutData{
                             {{0, 1, 2, 3, 4, 5, 6}}, {{6}}, {1, 2.5, 3, 5, 7, 8.1, 9}, 8.2, 10}));

INSTANTIATE_TEST_SUITE_P(
    TestInclusivityOfT0, TestSeismogramCut,
    testing::Values(TestSeismogramCutData{
        {{0, 1, 2, 3, 4, 5, 6}}, {{1, 2}}, {1, 2.5, 3, 5, 7, 8.1, 9}, 2.5, 3.5}));

INSTANTIATE_TEST_SUITE_P(
    TestKeepingInclusivityT1, TestSeismogramCut,
    testing::Values(TestSeismogramCutData{
        {{0, 1, 2, 3, 4, 5, 6}}, {{3, 4, 5}}, {1, 2.5, 3, 5, 7, 8.1, 9}, 4.5, 8.1}));