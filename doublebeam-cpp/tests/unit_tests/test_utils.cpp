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