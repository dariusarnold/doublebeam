#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include <complex>
#include <vector>

#include "fft.hpp"
#include "printing.hpp"

TEST(TestFFT, TestEmptyInput) {
    std::vector<double> in;
    FFT fft;
    std::vector<std::complex<double>> out;
    ASSERT_NO_THROW(out = fft.execute(in));
    ASSERT_EQ(out.size(), 0);
}


struct TestParam_t {
    std::vector<double> in;
    std::vector<std::complex<double>> expected;
};

class TestFFT : public testing::TestWithParam<TestParam_t> {
protected:
    FFT fft;
};

TEST_P(TestFFT, TestIfForwardTransformReturnsCorrectOutputByComparingWithNumpyResult) {
    // fftw doesn't return the redundant/symmetrical part, so the output array has size N/2 + 1
    // where N is the size of the input array.
    auto [in, expected] = GetParam();
    auto result = fft.execute(in);
    std::cout.precision(17);
    ASSERT_EQ(result.size(), expected.size()) << "Different size for expected and actual result.";
    for (auto i = 0; i < result.size(); ++i) {
        EXPECT_TRUE(AlmostEqual(result[i].real(), expected[i].real(), 14))
            << "Different real values at index " << i;
        EXPECT_TRUE(AlmostEqual(result[i].imag(), expected[i].imag(), 14))
            << "Different imaginary values at index " << i;
    }
}

using namespace std::complex_literals;
INSTANTIATE_TEST_SUITE_P(TestNormalValues, TestFFT,
                         testing::Values(TestParam_t{{0, 1, 2, 3},
                                                     {6. + 0.i, -2. + 2.i, -2. + 0.i}},
                                         TestParam_t{{0, 0, 0}, {0. + 0.i, 0. + 0.i}},
                                         TestParam_t{{2, 2}, {4. + 0.i, 0. + 0.i}},
                                         TestParam_t{{0, 1, 2, 1, 0, -1, -2, -1, 0},
                                                     {
                                                         0. + 0.i,
                                                         2.379385241571817 - 6.5373072234021485i,
                                                         -1.0320888862379565 + 1.2299956380506412i,
                                                         -1.5 + 0.8660254037844386i,
                                                         0.15270364466613895 - 0.02692577260715912i,
                                                     }}));
