#include <filesystem>

#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include "io.hpp"
#include "printing.hpp"
#include "seismodata.hpp"


std::ostream& operator<<(std::ostream& os, Seismogram<double> seismogram) {
    os << "Seismogram(t = ";
    for (auto t : seismogram.timesteps) {
        os << t << " ";
    }
    os << ", x = ";
    for (auto x : seismogram.data) {
        os << x << " ";
    }
    os << ")";
    return os;
}

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
    ASSERT_DOUBLE_EQ(s.timestep(), 0.004);
    // loading source and receiver file is already tested.
}

TEST_P(TestProjectLoading, TestRetrievingSeismograms) {
    SeismoData s(project_dir);
    // I replaced the first value of every seismogram with an increasing index to test if they are
    // loaded in the correct order.
    auto i = 1;
    for (const auto& source : s.sources()) {
        for (const auto& receiver : s.receivers()) {
            auto seismogram = s.get_seismogram(source, receiver);
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


struct TestCuttingData {
    TestCuttingData(double t0, double t1, size_t index0, size_t index1) :
            t0(t0), t1(t1), index0(index0), index1(index1) {}

    // cut times
    double t0, t1;
    // indices of expected range in seismogram, starting at 0. index1 is one past the last value,
    // so for a full cut of 1000 values specify index0=0, index1=1000
    size_t index0, index1;
};


class TestCutting : public testing::TestWithParam<TestCuttingData> {
protected:
    std::filesystem::path project_dir =
        current_source_path(__FILE__) / "data" / "sample_project_dir";
    SeismoData seismo_data{project_dir};

    // Read first seismogram which is used to test cutting
    std::vector<double> timesteps, amplitudes;
    TestCutting() :
            timesteps(read_timesteps(project_dir / "shotdata" / "source_001" / "receiver_001.txt")),
            amplitudes(
                read_amplitude(project_dir / "shotdata" / "source_001" / "receiver_001.txt")) {}

    // mock source/receiver so first seismogram from first source is extracted
    Source s{0, 0, 0, 1};
    Receiver r{0, 0, 0, 1};

    // helper function to pass source receiver so first seismogram is returned.
    // First seismogram is the test case for cutting, this saves writing s, r in every test.
    Seismogram<double> get_seismogram(double t0, double t1) {
        return seismo_data.get_seismogram(s, r, t0, t1);
    }
};

namespace testing::internal {
    bool operator==(const gsl::span<double> span, const std::vector<double>& vector) {
        if (span.size() != vector.size()) {
            return false;
        }
        for (auto i = 0; i < span.size(); ++i) {
            if (span[i] != vector[i]) {
                return false;
            }
        }
        return true;
    }
} // namespace testing::internal

TEST_P(TestCutting, TestKeepingFullSeismogram) {
    auto [t0, t1, index0, index1] = GetParam();
    auto expected_size = index1 - index0;
    auto out = get_seismogram(t0, t1);
    ASSERT_EQ(out.data.size(), out.timesteps.size()) << "Timesteps and data cut differently";
    ASSERT_EQ(out.size(), expected_size) << "Wrong size for cut seismogram.";
    gsl::span<double> expected_amplitudes(amplitudes.data() + index0, expected_size);
    gsl::span<double> expected_timesteps(timesteps.data() + index0, expected_size);
    EXPECT_EQ(out.data, expected_amplitudes) << "Wrong values for amplitude cut.";
    EXPECT_EQ(out.timesteps, expected_timesteps) << "Wrong values for timestep cut.";
}

INSTANTIATE_TEST_SUITE_P(TestEmptySeismogram, TestCutting,
                         testing::Values(TestCuttingData(0.0005, 0.000051, 0, 0)));

INSTANTIATE_TEST_SUITE_P(TestCuttingFullSeismogram, TestCutting,
                         testing::Values(TestCuttingData(-1, 5, 0, 1000)));

INSTANTIATE_TEST_SUITE_P(TestInclusivityOfT1, TestCutting,
                         testing::Values(TestCuttingData(0.002, 0.008, 1, 3)));

INSTANTIATE_TEST_SUITE_P(TestInclusivityOfT0, TestCutting,
                         testing::Values(TestCuttingData(0.004, 0.009, 1, 3)));

INSTANTIATE_TEST_SUITE_P(TestInclusivityOfT0AndT1, TestCutting,
                         testing::Values(TestCuttingData(0.004, 0.012, 1, 4)));

INSTANTIATE_TEST_SUITE_P(TestKeepingFirstValue, TestCutting,
                         testing::Values(TestCuttingData(-0.002, 0.002, 0, 1)));

INSTANTIATE_TEST_SUITE_P(TestKeepingLastValue, TestCutting,
                         testing::Values(TestCuttingData(3.997, 4.01, 999, 1000)));

INSTANTIATE_TEST_SUITE_P(TestCuttingGeneralCaseSmall, TestCutting,
                         testing::Values(TestCuttingData(0.008, 0.012, 2, 4)));

INSTANTIATE_TEST_SUITE_P(TestCuttingGeneralCaseMedium, TestCutting,
                         testing::Values(TestCuttingData(0.08, 0.12, 20, 31)));

INSTANTIATE_TEST_SUITE_P(TestCuttingGeneralCaseLarge, TestCutting,
                         testing::Values(TestCuttingData(0.8, 1.2, 200, 301)));