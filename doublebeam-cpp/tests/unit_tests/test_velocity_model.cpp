#include <gtest/gtest.h>

#include "model.h"


class TestInterfaceVelocity : public ::testing::Test {
protected:

    TestInterfaceVelocity() :
            vm(VelocityModel(std::vector<Layer>{
                    {0,   100, 1000, 0},
                    {100, 200, 1500, 1},
                    {200, 300, 3000, -1}})) {}

    VelocityModel vm;
};


TEST_F(TestInterfaceVelocity, TestEvaluation) {
    for (auto z: {99, 100, 101}) {
        auto[v_above, v_below] = vm.interface_velocities(z);
        EXPECT_DOUBLE_EQ(v_above, 1000);
        EXPECT_DOUBLE_EQ(v_below, 1600);
    }
    for (auto z: {199, 200, 201}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_DOUBLE_EQ(v_above, 1700);
        EXPECT_DOUBLE_EQ(v_below, 2800);
    }
}

TEST_F(TestInterfaceVelocity, TestVelocityTop) {
    /*
     * For depths above the layers midpoint, interface is between the top layer and air.
     * Velocity returned for air should be 0.
     */
    auto [v_above, v_below] = vm.interface_velocities(10.);
    EXPECT_DOUBLE_EQ(v_above, 0);
    EXPECT_DOUBLE_EQ(v_below, 1000);
}


TEST_F(TestInterfaceVelocity, TestVelocityBottom) {
    /**
     * For depths below the bottom layers mid point, return 0 as velocity
     * below bottom layer.
     */
    auto [v_above, v_below] = vm.interface_velocities(299.);
    EXPECT_DOUBLE_EQ(v_above, 2700);
    EXPECT_DOUBLE_EQ(v_below, 0);
}