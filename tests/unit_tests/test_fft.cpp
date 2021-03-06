/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include <complex>
#include <vector>

#include <gsl/span>
#include <gsl/span_ext>

#include "fft.hpp"
#include "printing.hpp"

TEST(TestFFT, TestEmptyInputVector) {
    std::vector<double> in;
    FFT fft;
    std::vector<std::complex<double>> out;
    ASSERT_NO_THROW(out = fft.execute(in));
    ASSERT_EQ(out.size(), 0);
}

TEST(TestFFT, TestEmptyInputRange) {
    std::vector<double> in;
    FFT fft;
    FFT::cvector_align out;
    ASSERT_NO_THROW(out = fft.execute(gsl::make_span(in)));
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

TEST_P(TestFFT, TestWithRanges) {
    auto [in, expected] = GetParam();
    auto result = fft.execute(gsl::make_span(in));
    std::cout.precision(17);
    ASSERT_EQ(result.size(), expected.size());
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
