#include <gtest/gtest.h>

#include "raytracing_helpers.hpp"


using depth_t = double;

class TestInterfaceCrossed : public testing::TestWithParam<depth_t> {
protected:
    InterfaceCrossed in_layer{10, 20};

};

class TestInInterface : public TestInterfaceCrossed {};
class TestOutsideInterface : public TestInterfaceCrossed {};

TEST_P(TestInInterface, TestDepthInsideLayer) {
    // test if depths inside the layer are recognized as having not left the layer
    auto depth = GetParam();
    auto state = state_type{0, 0, depth, 0, 0, 0};
    bool interface_crossed = in_layer(state);
    ASSERT_FALSE(interface_crossed);
}

TEST_P(TestOutsideInterface, TestDepthOutsideLayer) {
    // Test if depths outside the layer are recognized as having crossed the layer boundary
    auto depth = GetParam();
    auto state = state_type{0, 0, depth, 0, 0, 0};
    bool interface_crossed = in_layer(state);
    ASSERT_TRUE(interface_crossed);
}

INSTANTIATE_TEST_SUITE_P(InsideLayer, TestInInterface, testing::Values(10, 11, 15, 19, 20));

INSTANTIATE_TEST_SUITE_P(OutsideLayer, TestOutsideInterface, testing::Values(9.9999, 21, -5));