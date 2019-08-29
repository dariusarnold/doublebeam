#include "gtest/gtest.h"
#include <cmath>
#include "utils.h"

TEST(Radians, Equals) {
    EXPECT_DOUBLE_EQ(math::radians(0), 0.);
    EXPECT_DOUBLE_EQ(math::radians(360), 2*M_PI);
    EXPECT_DOUBLE_EQ(math::radians(360+180), 3*M_PI);
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

