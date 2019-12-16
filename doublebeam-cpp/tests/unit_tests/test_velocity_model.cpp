#include <utility>

#include <gtest/gtest.h>
#include <utils.hpp>

#include "model.hpp"


class TestInterfaceVelocity : public ::testing::Test {
protected:
    TestInterfaceVelocity() :
            vm(VelocityModel(std::vector<Layer>{{0, 100, 1000}, {100, 200, 1500}, {200, 300, 3000}},
                             1000, 1000)) {}

    VelocityModel vm;
};

TEST_F(TestInterfaceVelocity, TestEvaluation) {
    for (auto z : {99, 100, 101}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_DOUBLE_EQ(v_above, 1000) << "Error at depth " << z;
        EXPECT_DOUBLE_EQ(v_below, 1500) << "Error at depth " << z;
    }
    for (auto z : {199, 200, 201}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_DOUBLE_EQ(v_above, 1500) << "Error at depth " << z;
        EXPECT_DOUBLE_EQ(v_below, 3000) << "Error at depth " << z;
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
    EXPECT_DOUBLE_EQ(v_above, 3000);
    EXPECT_DOUBLE_EQ(v_below, 0);
}


// param is pair of depth, index. Index is the expected return value at depth.
class TestLayerIndex : public ::testing::TestWithParam<std::pair<double, size_t>> {
protected:
    TestLayerIndex() :
            vm(VelocityModel(std::vector<Layer>{{0, 100, 1800},
                                                {100, 200, 2400},
                                                {200, 300, 2400},
                                                {300, 400, 2700},
                                                {400, 500, 2250}},
                             1000, 1000)) {}

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
class TestLayerIndexReturnsInvalidOptional : public ::testing::Test {
protected:
    TestLayerIndexReturnsInvalidOptional() : vm(std::vector<Layer>{{0, 100, 0}}, 1000, 1000) {}

    VelocityModel vm;
};

TEST_F(TestLayerIndexReturnsInvalidOptional, TestDepthAboveModel) {
    ASSERT_FALSE(vm.layer_index(-1.1));
}

TEST_F(TestLayerIndexReturnsInvalidOptional, TestDepthBelowModel) {
    ASSERT_FALSE(vm.layer_index(1001));
}


TEST(TestCreateVelocityModelFromFile, TestSuccessfullRead) {
    // TODO better way to specify path, maybe mock file object
    std::filesystem::path filepath(
        "/home/darius/git/doublebeam/doublebeam-cpp/tests/unit_tests/data/model.txt");
    auto vm = read_velocity_file(filepath);
    VelocityModel expected({{0, 100, 1000}, {100, 200, 1200}}, 100, 100);
    EXPECT_EQ(vm, expected);
}

TEST(TestCreateVelocityModel, TestThrowOnInvalidWidths) {
    std::vector<Layer> l{{0, 1, 2}};
    EXPECT_THROW(VelocityModel(l, 100, 50), std::invalid_argument)
        << "Not throwing for model creation with different width along x and y axis when using "
           "default values for x0 and y0.";
    EXPECT_THROW(VelocityModel(l, 100, 100, 0, 1), std::invalid_argument)
        << "Not throwing for model creation with different width along x and y axis.";
}

TEST(TestCreateVelocityModel, TestThrowOnEmptyListOfLayers) {
    EXPECT_THROW(VelocityModel({}, 1, 1), std::invalid_argument)
        << "Not throwing when constructing from empty list of layers.";
}

TEST(TestCreateVelocityModel, DISABLED_TestThrowOnInvalidListOfLayers) {
    std::vector<Layer> l_gap{{0, 10, 1}, {20, 30, 1}};
    std::vector<Layer> l_overlapping{{0, 11, 1}, {10, 20, 1}};
    EXPECT_THROW(VelocityModel(l_gap, 1, 1), std::invalid_argument)
        << "Not throwing when constructing velocity model with gap between layers.";
    EXPECT_THROW(VelocityModel(l_overlapping, 1, 1), std::invalid_argument)
        << "Not throwing when constructing velocity model with overlapping layers.";
}


struct DepthVelocityData {
    double depth1;
    double depth2;
    double velocity;
};

std::ostream& operator<<(std::ostream& os, DepthVelocityData d) {
    os << "Depth1: " << d.depth1 << " m Depth2: " << d.depth2 << " m Velocity: " << d.velocity
       << "m/s";
    return os;
}

class TestTwoPointRayTracing_v_M : public ::testing::TestWithParam<DepthVelocityData> {
protected:
    // same model as Fang2019 fig. 3 but negative gradient in first layer
    TestTwoPointRayTracing_v_M() :
            model({{0, 100, 2600},
                   {100, 200, 2400},
                   {200, 300, 2400},
                   {300, 400, 2700},
                   {400, 500, 2250}},
                  1000, 1000) {}

    VelocityModel model;
};

TEST_P(TestTwoPointRayTracing_v_M, Test_v_M) {
    auto [depth1, depth2, velocity] = TestTwoPointRayTracing_v_M::GetParam();
    auto result = highest_velocity_between(depth1, depth2, model);
    EXPECT_EQ(result, velocity) << "Did not receive expected velocity " << velocity
                                << " between depth " << depth1 << " and " << depth2 << ".";
}

INSTANTIATE_TEST_SUITE_P(
    TestFindingHighestVelocityBetweenTwoDepths, TestTwoPointRayTracing_v_M,
    testing::Values(DepthVelocityData{0, 500, 2700}, DepthVelocityData{50, 500, 2700},
                    DepthVelocityData{50, 450, 2700}, DepthVelocityData{300, 400, 2700},
                    DepthVelocityData{50, 150, 2600}, DepthVelocityData{25, 150, 2600},
                    DepthVelocityData{301, 302, 2700}, DepthVelocityData{0, 1, 2600},
                    DepthVelocityData{50, 250, 2600}));

INSTANTIATE_TEST_SUITE_P(FixThisLater, TestTwoPointRayTracing_v_M,
                         testing::Values(DepthVelocityData{200, 300, 2400}));

struct DepthIntData {
    double z1, z2;
    size_t num_of_layers;
};

class TestNumberOfInterfaces : public ::testing::TestWithParam<DepthIntData> {
protected:
    // same model as Fang2019 fig. 3 but negative gradient in first layer
    TestNumberOfInterfaces() :
            model({{0, 100, 2600},
                   {100, 200, 2400},
                   {200, 300, 2400},
                   {300, 400, 2700},
                   {400, 500, 2250}},
                  1000, 1000) {}

    VelocityModel model;
};

std::vector<DepthIntData> values = {{0, 500, 4}, {50, 450, 4}, {50, 150, 1}, {200, 250, 0}};

INSTANTIATE_TEST_SUITE_P(TestInOrder, TestNumberOfInterfaces, testing::ValuesIn(values));


class TestModelGeometry : public testing::Test {
protected:
    TestModelGeometry() : model({{0, 500, 1}, {5000, 1000, 1}}, 1000, 1000) {}

    VelocityModel model;
};

TEST_F(TestModelGeometry, TestAboveModel) {
    ASSERT_FALSE(model.in_model(0, 0, -2));
}

TEST_F(TestModelGeometry, TestBelowModel) {
    ASSERT_FALSE(model.in_model(0, 0, 1001));
}

class TestModelGeometryParametrized
        : public TestModelGeometry,
          public testing::WithParamInterface<std::tuple<double, double, double>> {};

class TestInModel : public TestModelGeometryParametrized {};

TEST_P(TestInModel, TestIfPointInModel) {
    auto [x, y, z] = GetParam();
    EXPECT_TRUE(model.in_model(x, y, z)) << "Position " << (impl::Formatter(" ") << x << y << z)
                                         << " not recognized as inside model " << model;
}

class TestOutsideModel : public TestModelGeometryParametrized {};


TEST_P(TestOutsideModel, TestIfPointOutsideModel) {
    auto [x, y, z] = GetParam();
    EXPECT_FALSE(model.in_model(x, y, z)) << "Position " << (impl::Formatter(" ") << x << y << z)
                                          << "not recognized as outside model" << model;
}

INSTANTIATE_TEST_SUITE_P(TestCombinationsOfPointsInsideModel, TestInModel,
                         testing::Combine(testing::Values(0, 100, 900, 1000),
                                          testing::Values(0, 100, 900, 1000),
                                          testing::Values(0, 400, 500, 900, 1000)));

INSTANTIATE_TEST_SUITE_P(TestAboveModel, TestOutsideModel,
                         testing::Combine(testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(-1)));

INSTANTIATE_TEST_SUITE_P(TestBelowModel, TestOutsideModel,
                         testing::Combine(testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(1001)));

INSTANTIATE_TEST_SUITE_P(TestOutSideXLess, TestOutsideModel,
                         testing::Combine(testing::Values(-1),
                                          testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(0, 100, 500, 900, 1000)));

INSTANTIATE_TEST_SUITE_P(TestOutSideXGreater, TestOutsideModel,
                         testing::Combine(testing::Values(1001),
                                          testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(0, 100, 500, 900, 1000)));

INSTANTIATE_TEST_SUITE_P(TestOutSideYLess, TestOutsideModel,
                         testing::Combine(testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(-1),
                                          testing::Values(0, 100, 500, 900, 1000)));

INSTANTIATE_TEST_SUITE_P(TestOutSideYGreater, TestOutsideModel,
                         testing::Combine(testing::Values(0, 100, 500, 900, 1000),
                                          testing::Values(1001),
                                          testing::Values(0, 100, 500, 900, 1000)));