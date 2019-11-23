#include <gtest/gtest.h>

#include "raytracing_helpers.hpp"

TEST(TestDepthPredicate, TestInitialCallShouldBeFalse) {
    auto depth_pred = DepthPredicate(15);
    auto state = state_type{0, 0, 0, 0, 0, 0};
    bool result_of_first_call = depth_pred(state);
    ASSERT_FALSE(result_of_first_call);
}

TEST(TestDepthPredicate, TestIfWorking) {
    auto depth_pred = DepthPredicate(15);
    auto state = state_type{0, 0, 0, 0, 0, 0};
    // call once to set the current depth value
    depth_pred(state);
    state[Index::Z] = 10;
    bool result = depth_pred(state);
    EXPECT_FALSE(result) << "Predicate reported passing stop depth erroneously.";
    state[Index::Z] = 20;
    result = depth_pred(state);
    EXPECT_TRUE(result) << "Predicate did not recognize passing the stop depth.";
}
