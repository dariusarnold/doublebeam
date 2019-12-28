#include "units.hpp"
#include "utils.hpp"

#include "testing_utils.hpp"
#include "gtest/gtest.h"

#include <boost/container_hash/hash.hpp>

#include <cmath>
#include <unordered_set>


TEST(Radians, Equals) {
    EXPECT_DOUBLE_EQ(radians(0_deg).get(), 0.);
    EXPECT_DOUBLE_EQ(radians(360_deg).get(), 2 * M_PI);
    EXPECT_DOUBLE_EQ(radians(540_deg).get(), 3 * M_PI);
}


TEST(SameSign, EqualSign) {
    EXPECT_TRUE(math::same_sign(1.1, 2.9));
    EXPECT_TRUE(math::same_sign(9, 1));
    EXPECT_TRUE(math::same_sign(-4.2, -6.9));
    EXPECT_TRUE(math::same_sign(-10, -11));
}

TEST(SameSign, DifferentSigns) {
    EXPECT_FALSE(math::same_sign(-.1, 2.));
    EXPECT_FALSE(math::same_sign(2, -3));
}

TEST(SameSign, ExtremeNumbers) {
    EXPECT_TRUE(math::same_sign(1E302, 1E304));
    EXPECT_TRUE(math::same_sign(1E-300, 1E-323));
}

TEST(SameSign, ZeroTreatedAsPositive) {
    EXPECT_TRUE(math::same_sign(0, 1));
    EXPECT_TRUE(math::same_sign(1.1, 0.));
    EXPECT_FALSE(math::same_sign(-1, 0));
    EXPECT_FALSE(math::same_sign(0., -9.9));
}


struct AngleData {
    Radian angle;
    double x, y;
};

class TestClockWiseAngle : public ::testing::TestWithParam<AngleData> {
protected:
    double x_axis_x = 1;
    double x_axis_y = 0;
};

/**
 * The cartesian coordinate system:
 *  y
 *  ^
 *  |
 *  |
 *  |
 * -|------> x
 *  |
 */

TEST_P(TestClockWiseAngle, TestNormal) {
    // Test if function calculates angle between x axis and given vector correctly by
    // checking against manually computed test data.
    auto data = GetParam();
    auto angle_result = math::angle_clockwise(x_axis_x, x_axis_y, data.x, data.y);
    EXPECT_EQ(angle_result, data.angle.get()) << "For x: " << data.x << " y: " << data.y;
}

TEST_P(TestClockWiseAngle, TestSwapped) {
    // If the vectors are swapped, the "other angle" has to be returned. The value of this
    // other angle is 360° - first_angle.
    auto data = GetParam();
    auto angle_result = math::angle_clockwise(data.x, data.y, x_axis_x, x_axis_y);
    // module 360° since for 0° the function will not return 360°
    auto angle_expected = fmod(2 * M_PI - data.angle.get(), 2 * M_PI);
    EXPECT_EQ(angle_result, angle_expected) << "For x: " << data.x << " y: " << data.y;
}

TEST_F(TestClockWiseAngle, TestCase) {
    auto res = math::angle_clockwise(1, 0, -1, 1);
    EXPECT_EQ(degrees(Radian(res)).get(), (225_deg).get());
}

TEST_F(TestClockWiseAngle, TestCase2) {
    auto res = math::angle_clockwise(1, 0, 1, -1);
    EXPECT_EQ(degrees(Radian(res)).get(), (45_deg).get());
}

INSTANTIATE_TEST_SUITE_P(
    PREFIX, TestClockWiseAngle,
    ::testing::Values(AngleData{radians(0_deg), 1, 0}, AngleData{radians(315_deg), 1, 1},
                      AngleData{radians(270_deg), 0, 1}, AngleData{radians(225_deg), -1, 1},
                      AngleData{radians(180_deg), -1, 0}, AngleData{radians(135_deg), -1, -1},
                      AngleData{radians(90_deg), 0, -1}, AngleData{radians(45_deg), 1, -1}));


class TestRounding : public ::testing::Test {
protected:
    double value = 123.123456789123456789;
};

TEST_F(TestRounding, TestNormalCase) {
    EXPECT_EQ(123.12, math::round(value, 2));
    EXPECT_EQ(123.123, math::round(value, 3));
}

TEST_F(TestRounding, TestRoundingUp) {
    EXPECT_EQ(123.1235, math::round(value, 4));
}

TEST_F(TestRounding, TestCompleteLengthIsUnchanged) {
    EXPECT_EQ(value, math::round(value, 16));
}

TEST_F(TestRounding, TestZeroDigits) {
    EXPECT_EQ(static_cast<int>(value), math::round(value, 0));
}

class TestcumtrapzByComparingToScipy : public testing::Test {
protected:
    TestcumtrapzByComparingToScipy() {
        constexpr double start = -2, end = 2;
        constexpr static int steps = 20;
        constexpr double stepsize = std::abs(end - start) / (steps - 1);
        for (int i = 0; i < steps; ++i) {
            y.push_back(start + i * stepsize);
        }
    }
    // vector containing values of y = np.linspace(-2, 2, 20)
    std::vector<double> y;
};

TEST_F(TestcumtrapzByComparingToScipy, TestWithoutInitialParameterSet) {
    // compare results of my implementation with manually computed results from
    // scipy.integrate.cumtrapz(y)

    auto res = math::cumtrapz(y.begin(), y.end(), 0.2);
    std::vector<double> scipy_result{
        -3.7894736842105270e-01, -7.1578947368421053e-01, -1.0105263157894737e+00,
        -1.2631578947368420e+00, -1.4736842105263157e+00, -1.6421052631578947e+00,
        -1.7684210526315789e+00, -1.8526315789473684e+00, -1.8947368421052633e+00,
        -1.8947368421052633e+00, -1.8526315789473686e+00, -1.7684210526315791e+00,
        -1.6421052631578950e+00, -1.4736842105263159e+00, -1.2631578947368423e+00,
        -1.0105263157894737e+00, -7.1578947368421053e-01, -3.7894736842105270e-01,
        -5.5511151231257827e-17};
    ASSERT_EQ(res.size(), scipy_result.size()) << "Different size for actual and expected result.";
    for (auto i = 0; i < res.size(); ++i) {
        EXPECT_DOUBLE_EQ(res[i], scipy_result[i])
            << "Results differ at index " << i << ". expected " << res[i] << " got "
            << scipy_result[i];
    }
}

TEST_F(TestcumtrapzByComparingToScipy, TestInitialParameterSet) {
    // compare the results of my implementation with manually computed results for
    // cumtrapz(y, dx=0.2) from scipy.integrate.cumtrapz
    auto res = math::cumtrapz(y.begin(), y.end(), 0.2, 0.);
    std::vector<double> scipy_result{
        0.0000000000000000e+00,  -3.7894736842105270e-01, -7.1578947368421053e-01,
        -1.0105263157894737e+00, -1.2631578947368420e+00, -1.4736842105263157e+00,
        -1.6421052631578947e+00, -1.7684210526315789e+00, -1.8526315789473684e+00,
        -1.8947368421052633e+00, -1.8947368421052633e+00, -1.8526315789473686e+00,
        -1.7684210526315791e+00, -1.6421052631578950e+00, -1.4736842105263159e+00,
        -1.2631578947368423e+00, -1.0105263157894737e+00, -7.1578947368421053e-01,
        -3.7894736842105270e-01, -5.5511151231257827e-17};
    ASSERT_EQ(res.size(), scipy_result.size()) << "Different size for actual and expected result.";
    for (auto i = 0; i < res.size(); ++i) {
        EXPECT_DOUBLE_EQ(res[i], scipy_result[i])
            << "Results differ at index " << i << ". expected " << res[i] << " got "
            << scipy_result[i];
    }
}

TEST_F(TestcumtrapzByComparingToScipy, TestMultipleDistances) {
    // compare results from my implementation with manually computed results for x = y,
    // cumtrapz(y, x, initial=0)
    auto x = y;
    auto res = math::cumtrapz(y.begin(), y.end(), x.begin(), x.end(), 0.);
    std::vector<double> scipy_result{
        0.0000000000000000e+00,  -3.9889196675900279e-01, -7.5346260387811625e-01,
        -1.0637119113573408e+00, -1.3296398891966759e+00, -1.5512465373961219e+00,
        -1.7285318559556788e+00, -1.8614958448753465e+00, -1.9501385041551249e+00,
        -1.9944598337950141e+00, -1.9944598337950141e+00, -1.9501385041551249e+00,
        -1.8614958448753467e+00, -1.7285318559556793e+00, -1.5512465373961224e+00,
        -1.3296398891966761e+00, -1.0637119113573412e+00, -7.5346260387811725e-01,
        -3.9889196675900340e-01, -2.7755575615628914e-16};
    ASSERT_EQ(res.size(), scipy_result.size()) << "Different size for actual and expected result.";
    for (auto i = 0; i < res.size(); ++i) {
        EXPECT_DOUBLE_EQ(res[i], scipy_result[i])
            << "Results differ at index " << i << ". expected " << res[i] << " got "
            << scipy_result[i];
    }
}

TEST(Testcumtrapz, TestDynamicRayTracingUseCaseByComparingWithPythonResult) {
    std::vector<double> v_squared{3240000.0, 3240005.1194904763, 3240056.3145882203,
                                  3240568.284863763, 3245689.91708478};
    std::vector<double> travel_times{0.0, 2.101868808572426e-07, 2.312046556352558e-06,
                                     2.3329730060522415e-05, 0.00023341529691266796};
    std::vector<double> python_result{0.0, 0.6810060320023335, 7.491095943513512, 75.59495431382318,
                                      756.9295698666922};
    auto res = math::cumtrapz(v_squared.begin(), v_squared.end(), travel_times.begin(),
                              travel_times.end(), 0);
    ASSERT_EQ(res.size(), python_result.size()) << "Different size for actual and expected result.";
    for (auto i = 0; i < res.size(); ++i) {
        EXPECT_DOUBLE_EQ(res[i], python_result[i])
            << "Results differ at index " << i << ". expected " << res[i] << " got "
            << python_result[i];
    }
}

TEST(TestAngle, Test90Degrees) {
    EXPECT_EQ(math::angle(1, 0, 0, 0, 1, 0), radians(90_deg));
    EXPECT_EQ(math::angle(1, 0, 0, 0, 1, 0, false), radians(90_deg));
}

TEST(TestAngle, TestParallel) {
    EXPECT_EQ(math::angle(.8, .2, .4, .8, .2, .4), Radian(0));
    EXPECT_EQ(math::angle(.8, .2, .4, .8, .2, .4, false), radians(180_deg));
}

TEST(TestAngle, TestFloatClippingForAcosLargerOne) {
    EXPECT_EQ(math::angle(0.874469283050132, 0.262553720250597, 0.477968688795641,
                          0.874469283050132, 0.262553720250597, 0.477968688795641),
              Radian(0));
    EXPECT_EQ(math::angle(0.874469283050132, 0.262553720250597, 0.477968688795641,
                          0.874469283050132, 0.262553720250597, 0.477968688795641, false),
              radians(180_deg));
}

TEST(TestAngle, OrderOfInputVectorsShouldntMatter) {
    EXPECT_EQ(math::angle(2, 1, 0, 1, 0, 0), math::angle(1, 0, 0, 2, 1, 0));
    EXPECT_EQ(math::angle(2, 1, 0, 1, 0, 0, false), math::angle(1, 0, 0, 2, 1, 0, false));
}

TEST(TestIndicesFromRayCode, TestCase) {
    std::vector<int> expected_indices{0, 1, 2, 3, 4, 4, 3, 3, 3, 2, 1, 0};
    auto indices = seismo::ray_code_to_layer_indices("TTTTRTRRTTT", 1);
    ASSERT_EQ(indices.size(), expected_indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        EXPECT_EQ(indices[i], expected_indices[i]) << "Difference at index " << i;
    }
}
struct LinspaceExpectedResult {
    double start, stop;
    int num;
    std::vector<double> res;
};

std::ostream& operator<<(std::ostream& os, const LinspaceExpectedResult& l) {
    os << "linspace(start=" << l.start << ", stop=" << l.stop << ", num=" << l.num << ") = ";
    std::copy(l.res.begin(), l.res.end(), std::ostream_iterator<double>(os, " "));
    return os;
}

class TestLinspaceByComparingWithExpectedResult
        : public testing::TestWithParam<LinspaceExpectedResult> {};

TEST_P(TestLinspaceByComparingWithExpectedResult, CompareResults) {
    auto params = GetParam();
    auto a = math::linspace(params.start, params.stop, params.num);
    ASSERT_EQ(a.size(), params.res.size()) << "Different size for result.";
    ASSERT_EQ(a, params.res) << "Different array.";
}

INSTANTIATE_TEST_SUITE_P(
    TestLinspaceIncreasing, TestLinspaceByComparingWithExpectedResult,
    testing::Values(LinspaceExpectedResult{0, 10, 0, {}}, LinspaceExpectedResult{0, 1, 1, {0}},
                    LinspaceExpectedResult{0, 10, 2, {0, 10}},
                    LinspaceExpectedResult{0, 10, 11, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}},
                    LinspaceExpectedResult{1, 2, 6, {1, 1.2, 1.4, 1.6, 1.8, 2.0}},
                    LinspaceExpectedResult{-2, 2, 5, {-2, -1, 0, 1, 2}}));

INSTANTIATE_TEST_SUITE_P(
    TestLinspaceDecreasing, TestLinspaceByComparingWithExpectedResult,
    testing::Values(LinspaceExpectedResult{10, 0, 0, {}}, LinspaceExpectedResult{1, 0, 1, {1}},
                    LinspaceExpectedResult{10, 0, 2, {10, 0}},
                    LinspaceExpectedResult{10, 0, 11, {10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}},
                    LinspaceExpectedResult{2, 1, 6, {2.0, 1.8, 1.6, 1.4, 1.2, 1.0}},
                    LinspaceExpectedResult{2, -2, 5, {2, 1, 0, -1, -2}}));


TEST(TestFormatter, TestWithoutSeparator) {
    auto s = "Welt";
    EXPECT_EQ(std::string(impl::Formatter() << "Hallo " << s << 123), "Hallo Welt123");
}

TEST(TestFormatter, TestWithSeparator) {
    EXPECT_EQ(std::string(impl::Formatter(", ") << 1.9 << 2 << 3.0), "1.9, 2, 3");
}

TEST(TestFormatter, TestToOStream) {
    std::stringstream os;
    os << (impl::Formatter() << "Hallo Welt");
    EXPECT_EQ(os.str(), "Hallo Welt");
}

TEST(TestSlowness3D, TestIfVelocityRelationIsFullfilled) {
    // slowness is defined as p = 1/v, so velocity should be recovered by v = 1/p
    auto [px, py, pz] =
        seismo::slowness_3D(radians(20_deg), radians(10_deg), 1000_meter_per_second);
    EXPECT_DOUBLE_EQ(1000., 1 / math::length(px.get(), py.get(), pz.get()));
}

TEST(TestSlowness3D, TestUpgoingSlownessShouldBeNegative) {
    auto [px, py, pz] = seismo::slowness_3D(radians(170_deg), 0_rad, 1000_meter_per_second);
    ASSERT_LT(pz.get(), 0);
}

TEST(TestSlowness3D, HorizontalRayShouldHaveZeroAsVerticalSlowness) {
    auto [px, py, pz] =
        seismo::slowness_3D(radians(90_deg), radians(10_deg), 1000_meter_per_second);
    EXPECT_TRUE(AlmostEqual(pz.get(), 0., 17));
    EXPECT_NE(px.get(), py.get());
    EXPECT_NE(px.get(), 0);
    EXPECT_NE(py.get(), 0);
}

TEST(TestGridCoordinates, TestEmptySizes) {
    auto result = seismo::grid_coordinates(0_meter, 1_meter, 0_meter, 1_meter, 1_meter, 0, 0);
    ASSERT_EQ(result.size(), 0);
}

TEST(TestGridCoordinates, TestSizeOneGrid) {
    auto result = seismo::grid_coordinates(0_meter, 1_meter, 0_meter, 1_meter, 1_meter, 1, 1);
    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result[0], Position(0_meter, 0_meter, 1_meter));
}

TEST(TestGridCoordinates, TestNormalUseCase) {
    auto result = seismo::grid_coordinates(0_meter, 3_meter, 4_meter, 6_meter, 1_meter, 4, 3);
    ASSERT_EQ(result.size(), 4 * 3);
    std::unordered_set<Position, boost::hash<Position>> coordinates;
    for (auto x : {0_meter, 1_meter, 2_meter, 3_meter}) {
        for (auto y : {4_meter, 5_meter, 6_meter}) {
            coordinates.emplace(x, y, 1_meter);
        }
    }
    for (auto pos : result) {
        EXPECT_EQ(coordinates.count(pos), 1);
    }
}

TEST(TestGridCoordinates, TestNormalUseCaseWithReversedGridExtent) {
    auto result = seismo::grid_coordinates(3_meter, 0_meter, 6_meter, 4_meter, 1_meter, 4, 3);
    ASSERT_EQ(result.size(), 4 * 3);
    std::unordered_set<Position, boost::hash<Position>> coordinates;
    for (auto x : {0_meter, 1_meter, 2_meter, 3_meter}) {
        for (auto y : {4_meter, 5_meter, 6_meter}) {
            coordinates.emplace(x, y, 1_meter);
        }
    }
    for (auto pos : result) {
        EXPECT_EQ(coordinates.count(pos), 1);
    }
}


TEST(TestGoertzel, TestFloatNumbers) {
    std::vector<float> input{0, 1, 2, 3, 2, 1, 0};
    using namespace std::complex_literals;
    std::vector<std::complex<float>> expected_result{9.f + 0.if,
                                                     -4.5489173f - 2.190643if,
                                                     0.19202155f + 0.24078733if,
                                                     -0.14310412f - 0.6269801if,
                                                     -0.14310412f + 0.6269801if,
                                                     0.19202155f - 0.24078733if,
                                                     -4.5489173f + 2.190643if};
    for (auto bin = 0; bin < input.size(); ++bin) {
        EXPECT_TRUE(Close(math::goertzel(input, bin), expected_result[bin], {1E-5, 1E-5}))
            << "Different result in bin " << bin << ".";
    }
}

TEST(TestGoertzel, TestEmptyDataThrows) {
    std::vector<double> empty;
    ASSERT_THROW(math::goertzel(empty, 0), std::invalid_argument);
}

TEST(TestGoertzel, TestInvalidFrequencyBinThrows) {
    ASSERT_THROW(math::goertzel(std::vector<double>{1, 2, 3}, 5), std::invalid_argument);
}

TEST(TestGoertzel, TestDoubleNumbers) {
    std::vector<double> input{0, 1, 2, 3, 2, 1, 0};
    using namespace std::complex_literals;
    std::vector<std::complex<double>> expected_result{9. + 0.i,
                                                      -4.548917339522305 - 2.1906431337674115i,
                                                      0.19202147163009653 + 0.24078730940376436i,
                                                      -0.14310413210778994 - 0.6269801688313521i,
                                                      -0.14310413210778994 + 0.6269801688313521i,
                                                      0.19202147163009653 - 0.24078730940376436i,
                                                      -4.548917339522305 + 2.1906431337674115i};
    for (auto bin = 0; bin < input.size(); ++bin) {
        EXPECT_TRUE(Close(math::goertzel(input, bin), expected_result[bin]))
            << "Different result in bin " << bin << ".";
    }
}

TEST(TestGoertzel, TestClosestFrequency) {
    std::vector<double> input{0, 1, 2, 3, 2, 1, 0};
    using namespace std::complex_literals;
    //
    std::vector<std::complex<double>> expected_result{9. + 0i,
                                                      -4.5489173395223084 + -2.1906431337674093i,
                                                      0.19202147163009603 + 0.24078730940376419i,
                                                      -0.14310413210779105 + -0.62698016883135144i,
                                                      -0.1431041321077906 + 0.62698016883135166i,
                                                      0.19202147163009556 + -0.24078730940376414i,
                                                      -4.5489173395223013,
                                                      2.1906431337674155i};
    double sample_rate_s = 0.04;
    AngularFrequency sampling_frequency(2 * M_PI / sample_rate_s);
    // assuming input is sampled with sampling rate of 0.04 seconds this would give us
    // 0., 3.57142857, 7.14285714, 10.71428571 as the frequency bins (Hz) of the result.
    EXPECT_EQ(math::fft_closest_frequency(input, 0._rad_per_sec, sampling_frequency),
              expected_result[0]);
    EXPECT_EQ(
        math::fft_closest_frequency(input, AngularFrequency(2 * M_PI * 1.5), sampling_frequency),
        expected_result[0]);
    EXPECT_EQ(
        math::fft_closest_frequency(input, AngularFrequency(2 * M_PI * 2.5), sampling_frequency),
        expected_result[1]);
    EXPECT_EQ(
        math::fft_closest_frequency(input, AngularFrequency(2 * M_PI * 5.5), sampling_frequency),
        expected_result[2]);
    EXPECT_EQ(
        math::fft_closest_frequency(input, AngularFrequency(2 * M_PI * 11.), sampling_frequency),
        expected_result[3]);
}

TEST(TestGoertzel, TestSingleRealInputValue) {
    std::vector<double> input{42};
    std::complex<double> expected_result{42, 0};
    auto result = math::goertzel(input, 0);
    EXPECT_EQ(result, expected_result);
}

TEST(TestGoertzel, TestTwoRealInputValue) {
    std::vector<double> input{42, 24};
    std::vector<std::complex<double>> expected_result{{66, 0}, {18, 0}};
    for (auto bin = 0; bin < input.size(); ++bin) {
        EXPECT_TRUE(Close(math::goertzel(input, bin), expected_result[bin]))
            << "Different result in bin " << bin << ".";
    }
}


TEST(TestBetween, TestInclusivityLeft) {
    bool result = math::between(0, 0, 1);
    EXPECT_TRUE(result);
}

TEST(TestBetween, TestInclusivityRight) {
    bool result = math::between(0, 1, 1);
    EXPECT_TRUE(result);
}

TEST(TestBetween, TestNormalValueInside) {
    bool result = math::between(0., 0.96, 1.);
    EXPECT_TRUE(result);
}

TEST(TestBetween, TestNormalValueOutsideRight) {
    bool result = math::between(0., 1.01, 1.);
    EXPECT_FALSE(result);
}

TEST(TestBetween, TestNormalValueOutsideLeft) {
    bool result = math::between(0., -0.01, 1.);
    EXPECT_FALSE(result);
}

TEST(TestBetween, TestNormalInclusivityBothSides) {
    bool result = math::between(4.2, 4.2, 4.2);
    EXPECT_TRUE(result);
}

TEST(TestMatrixInverse, TestCase1) {
    auto [i11, i12, i13, i21, i22, i23, i31, i32, i33] =
        math::inv(1., 2., 0., 3., 0., 7., 0., 6., 5.);
    auto [e11, e12, e13, e21, e22, e23, e31, e32, e33] =
        std::make_tuple(0.5833333333333333333, 0.13888888888888888889, -0.19444444444444444443,
                        0.20833333333333333335, -0.069444444444444444444, 0.097222222222222222216,
                        -0.25, 0.083333333333333333333, 0.083333333333333333333);
    // I think this can be done better but I dont know how since get expects a constant expression
    // index.
    EXPECT_FLOAT_EQ(i11, e11);
    EXPECT_FLOAT_EQ(i12, e12);
    EXPECT_FLOAT_EQ(i13, e13);
    EXPECT_FLOAT_EQ(i21, e21);
    EXPECT_FLOAT_EQ(i22, e22);
    EXPECT_FLOAT_EQ(i23, e23);
    EXPECT_FLOAT_EQ(i31, e31);
    EXPECT_FLOAT_EQ(i32, e32);
    EXPECT_FLOAT_EQ(i33, e33);
}

TEST(TestVectorScaling, TestMakingVectorLonger) {
    auto [x, y, z] = math::scale_vector({1, 0, 0}, 2.5);
    EXPECT_DOUBLE_EQ(x, 2.5);
    EXPECT_EQ(y, 0);
    EXPECT_EQ(z, 0);
}

TEST(TestVectorScaling, TestMakingVectorShorter) {
    auto [x, y, z] = math::scale_vector({1, 0, 0}, .55);
    EXPECT_DOUBLE_EQ(x, .55);
    EXPECT_EQ(y, 0);
    EXPECT_EQ(z, 0);
}

struct TestCrossProdcutData {
    std::tuple<double, double, double> vector1;
    std::tuple<double, double, double> vector2;
    std::tuple<double, double, double> result_vector;
};

class TestCrossProduct : public testing::TestWithParam<TestCrossProdcutData> {};

TEST_P(TestCrossProduct, TestByComparingWithGivenValidResult) {
    auto data = GetParam();
    auto [x1, y1, z1] = data.vector1;
    auto [x2, y2, z2] = data.vector2;
    auto [result_x, result_y, result_z] = data.result_vector;
    auto [x, y, z] = math::cross(x1, y1, z1, x2, y2, z2);

    EXPECT_EQ(x, result_x);
    EXPECT_EQ(y, result_y);
    EXPECT_EQ(z, result_z);
}

INSTANTIATE_TEST_SUITE_P(TestUnitVectors, TestCrossProduct,
                         testing::Values(TestCrossProdcutData{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                         TestCrossProdcutData{{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},
                                         TestCrossProdcutData{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}));


TEST(CreateVectorArc, TestSimpleThreeVectorCase) {
    math::Vector2 central_direction{1, 1};
    auto number_of_vectors_to_generate = 3;
    auto vectors = math::generate_vector_arc(number_of_vectors_to_generate, central_direction);
    ASSERT_EQ(vectors.size(), number_of_vectors_to_generate);
    // test left vector
    EXPECT_DOUBLE_EQ(
        math::angle_clockwise(central_direction.x, central_direction.y, vectors[0].x, vectors[0].y),
        radians(270_deg).get());
    // test middle vector
    EXPECT_TRUE(Close(
        math::angle(central_direction.x, central_direction.y, 0, vectors[1].x, vectors[1].y, 0)
            .get(),
        radians(0_deg).get(), 0., 1.5e-8));
    // test right vector
    EXPECT_DOUBLE_EQ(
        math::angle_clockwise(central_direction.x, central_direction.y, vectors[2].x, vectors[2].y),
        radians(90_deg).get());
}

TEST(CreateVectorArc, TestCaseWithMoreVectors) {
    math::Vector2 central_direction{-1, 1};
    auto number_of_vectors_to_generate = 12;
    auto vectors = math::generate_vector_arc(number_of_vectors_to_generate, central_direction);
    ASSERT_EQ(vectors.size(), number_of_vectors_to_generate);
    EXPECT_DOUBLE_EQ(
        math::angle_clockwise(central_direction.x, central_direction.y, vectors[0].x, vectors[0].y),
        radians(270_deg).get())
        << " Leftmost vector is not 90 degrees to central direction.";
    EXPECT_DOUBLE_EQ(math::angle_clockwise(central_direction.x, central_direction.y,
                                           vectors.back().x, vectors.back().y),
                     radians(90_deg).get())
        << "Rightmost vector is not 90 degrees to central direction.";
    Degree expected_angle(180. / (number_of_vectors_to_generate - 1));
    for (auto i = 1; i < vectors.size(); ++i) {
        Radian angle(
            math::angle(vectors[i - 1].x, vectors[i - 1].y, 0, vectors[i].x, vectors[i].y, 0));
        EXPECT_TRUE(Close(angle.get(), radians(expected_angle).get(), 0., 1e-10))
            << "Vectors with index " << i - 1 << " and " << i << " have angle of " << degrees(angle)
            << ", but it should be " << expected_angle;
    }
}