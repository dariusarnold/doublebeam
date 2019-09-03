#include "model.hpp"
#include "twopoint.hpp"
#include "gtest/gtest.h"


class TestTwopointRayTracing : public ::testing::Test {
protected:
    TestTwopointRayTracing() : twopoint(model) {}

    VelocityModel model = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    TwoPointRayTracing twopoint;
};

TEST_F(TestTwopointRayTracing, TestThrowWhenOutOfRange) {
    auto [top, bottom] = model.get_top_bottom();

    EXPECT_THROW(twopoint.trace({0, 0, top - 1}, {0, 0, bottom}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, bottom + 1}, {0, 0, bottom}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, top}, {0, 0, top - 1}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, top}, {0, 0, bottom + 1}), std::domain_error);
}