#include <utility>

#include <gtest/gtest.h>

#include "model.hpp"


class TestInterfaceVelocity : public ::testing::Test {
protected:
    TestInterfaceVelocity()
    : vm(VelocityModel(
          std::vector<Layer>{{0, 100, 1000, 0}, {100, 200, 1500, 1}, {200, 300, 3000, -1}})) {}

    VelocityModel vm;
};

TEST_F(TestInterfaceVelocity, TestEvaluation) {
    for (auto z : {99, 100, 101}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_DOUBLE_EQ(v_above, 1000);
        EXPECT_DOUBLE_EQ(v_below, 1600);
    }
    for (auto z : {199, 200, 201}) {
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


class TestInModel : public TestInterfaceVelocity {
protected:
};

TEST_F(TestInModel, TestAboveModel) {
    EXPECT_FALSE(vm.in_model(-0.1));
}

TEST_F(TestInModel, TestBelowModel) {
    EXPECT_FALSE(vm.in_model(300.1));
}

TEST_F(TestInModel, TestInterfacesShouldBeInModel) {
    for (double interface_depth : {0, 100, 200, 300}) {
        EXPECT_TRUE(vm.in_model(interface_depth)) << "depth " << interface_depth;
    }
}


// param is pair of depth, index. Index is the expected return value at depth.
class TestLayerIndex : public ::testing::TestWithParam<std::pair<double, size_t>> {
protected:
    TestLayerIndex()
    : vm(VelocityModel(std::vector<Layer>{{0, 100, 1800, 4},
                                          {100, 200, 2400, 0},
                                          {200, 300, 2400, 1},
                                          {300, 400, 2700, 0},
                                          {400, 500, 2250, 1.5}})) {}

    VelocityModel vm;
};

TEST_P(TestLayerIndex, TestCorrectIndex) {
    auto [depth, index] = GetParam();
    EXPECT_EQ(vm.layer_index(depth), index);
}

/*
 * Test if indices in the middle of the layer are returned correctly.
 */
INSTANTIATE_TEST_SUITE_P(TestIndexInLayers, TestLayerIndex,
                         testing::Values(std::make_pair(50, 0), std::make_pair(140, 1),
                                         std::make_pair(250, 2), std::make_pair(350, 3),
                                         std::make_pair(450, 4)));
/*
 * Upper border of layer is inclusive, lower border is exclusive.
 * Bottom layer is special case where lower border in inclusive as well.
 */
INSTANTIATE_TEST_SUITE_P(TestIndexOnInterface, TestLayerIndex,
                         testing::Values(std::make_pair(0, 0), std::make_pair(100, 1),
                                         std::make_pair(400, 4), std::make_pair(500, 4)));

/*
 * Test if access to layer index outside of model boundaries throws error
 */
class TestLayerIndexThrows : public ::testing::Test {
protected:
    TestLayerIndexThrows() : vm(std::vector<Layer>{{0, 100, 0, 0}}) {}

    VelocityModel vm;
};

TEST_F(TestLayerIndexThrows, TestDepthAboveModel) {
    ASSERT_THROW(vm.layer_index(-1.1), std::domain_error);
}

TEST_F(TestLayerIndexThrows, TestDepthBelowModel) {
    ASSERT_THROW(vm.layer_index(1001), std::domain_error);
}


// inherit to get the same velocity model but different test names
class TestInterfaceVelocities : public TestLayerIndex {};

TEST_F(TestInterfaceVelocities, NoInterfaceBetweenDepths) {
    // interval should be empty
    auto [begin, end] = vm.interface_velocities(50, 50);
    EXPECT_EQ(begin, end);
}

TEST_F(TestInterfaceVelocities, OneInterfaceBetweenDepths) {
    // interval should contain the two values above/below the enclosed interface
    auto [begin, end] = vm.interface_velocities(50, 150);
    EXPECT_EQ(std::distance(begin, end), 2);
    EXPECT_EQ(*begin, 2200);
    EXPECT_EQ(*(begin + 1), 2400);
}

TEST_F(TestInterfaceVelocities, TwoInterfacesBetweenDepths) {
    // interval should contain four values from two above/below combinations of the two interfaces
    auto [begin, end] = vm.interface_velocities(150, 350);
    EXPECT_EQ(std::distance(begin, end), 4);
    for (auto velocity : {2400, 2600, 2700, 2700}) {
        EXPECT_EQ(*begin++, velocity);
    }
}


TEST(TestCreateVelocityModelFromFile, TestSuccessfullRead) {
    // TODO better way to specify path, maybe mock file object
    std::filesystem::path filepath(
        "/home/darius/git/doublebeam/doublebeam-cpp/tests/unit_tests/data/model.txt");
    auto vm = read_velocity_file(filepath);
    VelocityModel expected({{0, 100, 1000, 1}, {100, 200, 1200, -1}});
    EXPECT_EQ(vm, expected);
}