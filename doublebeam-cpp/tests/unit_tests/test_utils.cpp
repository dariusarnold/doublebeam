#include "utils.hpp"
#include "gtest/gtest.h"
#include <cmath>


TEST(Radians, Equals) {
    EXPECT_DOUBLE_EQ(math::radians(0), 0.);
    EXPECT_DOUBLE_EQ(math::radians(360), 2 * M_PI);
    EXPECT_DOUBLE_EQ(math::radians(360 + 180), 3 * M_PI);
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
    double angle;
    double x, y;
};

class TestClockWiseAngle : public ::testing::TestWithParam<AngleData> {
protected:
    double x_axis_x = 1;
    double x_axis_y = 0;
};

TEST_P(TestClockWiseAngle, TestNormal) {
    // Test if function calculates angle between x axis and given vector correctly by
    // checking against manually computed test data.
    auto data = GetParam();
    auto angle_result = math::angle_clockwise(x_axis_x, x_axis_y, data.x, data.y);
    EXPECT_EQ(angle_result, data.angle) << "For x: " << data.x << " y: " << data.y;
}

TEST_P(TestClockWiseAngle, TestSwapped) {
    // If the vectors are swapped, the "other angle" has to be returned. The value of this
    // other angle is 360° - first_angle.
    auto data = GetParam();
    auto angle_result = math::angle_clockwise(data.x, data.y, x_axis_x, x_axis_y);
    // module 360° since for 0° the function will not return 360°
    auto angle_expected = fmod(2 * M_PI - data.angle, 2 * M_PI);
    EXPECT_EQ(angle_result, angle_expected) << "For x: " << data.x << " y: " << data.y;
}

INSTANTIATE_TEST_SUITE_P(
    PREFIX, TestClockWiseAngle,
    ::testing::Values(AngleData{math::radians(0), 1, 0}, AngleData{math::radians(315), 1, 1},
                      AngleData{math::radians(270), 0, 1}, AngleData{math::radians(225), -1, 1},
                      AngleData{math::radians(180), -1, 0}, AngleData{math::radians(135), -1, -1},
                      AngleData{math::radians(90), 0, -1}, AngleData{math::radians(45), 1, -1}));


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