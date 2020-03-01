/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include <utility>

#include <gtest/gtest.h>
#include <utils.hpp>

#include "model.hpp"


class TestInterfaceVelocity : public ::testing::Test {
protected:
    TestInterfaceVelocity() :
            vm(VelocityModel(std::vector<Layer>{{0_meter, 100_meter, 1000_meter_per_second},
                                                {100_meter, 200_meter, 1500_meter_per_second},
                                                {200_meter, 300_meter, 3000_meter_per_second}},
                             1000_meter, 1000_meter)) {}

    VelocityModel vm;
};

TEST_F(TestInterfaceVelocity, TestEvaluation) {
    for (auto z : {99_meter, 100_meter, 101_meter}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_EQ(v_above, 1000_meter_per_second) << "Error at depth " << z;
        EXPECT_EQ(v_below, 1500_meter_per_second) << "Error at depth " << z;
    }
    for (auto z : {199_meter, 200_meter, 201_meter}) {
        auto [v_above, v_below] = vm.interface_velocities(z);
        EXPECT_EQ(v_above, 1500_meter_per_second) << "Error at depth " << z;
        EXPECT_EQ(v_below, 3000_meter_per_second) << "Error at depth " << z;
    }
}

TEST_F(TestInterfaceVelocity, TestVelocityTop) {
    /*
     * For depths above the layers midpoint, interface is between the top layer and air.
     * Velocity returned for air should be 0.
     */
    auto [v_above, v_below] = vm.interface_velocities(10._meter);
    EXPECT_EQ(v_above, 0_meter_per_second);
    EXPECT_EQ(v_below, 1000_meter_per_second);
}

TEST_F(TestInterfaceVelocity, TestVelocityBottom) {
    /**
     * For depths below the bottom layers mid point, return 0 as velocity
     * below bottom layer.
     */
    auto [v_above, v_below] = vm.interface_velocities(299._meter);
    EXPECT_EQ(v_above, 3000_meter_per_second);
    EXPECT_EQ(v_below, 0_meter_per_second);
}


// param is pair of depth, index. Index is the expected return value at depth.
class TestLayerIndex : public ::testing::TestWithParam<std::pair<Meter, size_t>> {
protected:
    TestLayerIndex() :
            vm(VelocityModel(std::vector<Layer>{{0_meter, 100_meter, 1800_meter_per_second},
                                                {100_meter, 200_meter, 2400_meter_per_second},
                                                {200_meter, 300_meter, 2400_meter_per_second},
                                                {300_meter, 400_meter, 2700_meter_per_second},
                                                {400_meter, 500_meter, 2250_meter_per_second}},
                             1000_meter, 1000_meter)) {}

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
                         testing::Values(std::make_pair(50_meter, 0), std::make_pair(140_meter, 1),
                                         std::make_pair(250_meter, 2), std::make_pair(350_meter, 3),
                                         std::make_pair(450_meter, 4)));
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
    TestLayerIndexReturnsInvalidOptional() :
            vm(std::vector<Layer>{{0_meter, 100_meter, 0_meter_per_second}}, 1000_meter,
               1000_meter) {}

    VelocityModel vm;
};

TEST_F(TestLayerIndexReturnsInvalidOptional, TestDepthAboveModel) {
    ASSERT_FALSE(vm.layer_index(-1.1_meter));
}

TEST_F(TestLayerIndexReturnsInvalidOptional, TestDepthBelowModel) {
    ASSERT_FALSE(vm.layer_index(1001_meter));
}


TEST(TestCreateVelocityModelFromFile, TestSuccessfullRead) {
    // TODO better way to specify path, maybe mock file object
    std::filesystem::path filepath(
        "/home/darius/git/doublebeam/tests/unit_tests/data/model.txt");
    auto vm = read_velocity_file(filepath);
    VelocityModel expected({{0_meter, 100_meter, 1000_meter_per_second},
                            {100_meter, 200_meter, 1200_meter_per_second}},
                           100_meter, 100_meter);
    EXPECT_EQ(vm, expected);
}

TEST(TestCreateVelocityModel, TestThrowOnEmptyListOfLayers) {
    EXPECT_THROW(VelocityModel({}, 1_meter, 1_meter), std::invalid_argument)
        << "Not throwing when constructing from empty list of layers.";
}

TEST(TestCreateVelocityModel, DISABLED_TestThrowOnInvalidListOfLayers) {
    std::vector<Layer> l_gap{{0_meter, 10_meter, 1_meter_per_second},
                             {20_meter, 30_meter, 1_meter_per_second}};
    std::vector<Layer> l_overlapping{{0_meter, 11_meter, 1_meter_per_second},
                                     {10_meter, 20_meter, 1_meter_per_second}};
    EXPECT_THROW(VelocityModel(l_gap, 1_meter, 1_meter), std::invalid_argument)
        << "Not throwing when constructing velocity model with gap between layers.";
    EXPECT_THROW(VelocityModel(l_overlapping, 1_meter, 1_meter), std::invalid_argument)
        << "Not throwing when constructing velocity model with overlapping layers.";
}


struct DepthVelocityData {
    DepthVelocityData(double d1, double d2, double vel) : depth1(d1), depth2(d2), velocity(vel) {}
    Meter depth1;
    Meter depth2;
    Velocity velocity;
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
            model({{0_meter, 100_meter, 2600_meter_per_second},
                   {100_meter, 200_meter, 2400_meter_per_second},
                   {200_meter, 300_meter, 2400_meter_per_second},
                   {300_meter, 400_meter, 2700_meter_per_second},
                   {400_meter, 500_meter, 2250_meter_per_second}},
                  1000_meter, 1000_meter) {}

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

// In the case of this function, the bottom of a layer should not count as belonging to the
// layer below. Instead when asking for the velocity between top and bottom depth of a layer,
// only the layer should be considered.
INSTANTIATE_TEST_SUITE_P(DISABLED_FixThisLater, TestTwoPointRayTracing_v_M,
                         testing::Values(DepthVelocityData{200, 300, 2400}));

struct DepthIntData {
    Meter z1, z2;
    size_t num_of_layers;
};

class TestNumberOfInterfaces : public ::testing::TestWithParam<DepthIntData> {
protected:
    // same model as Fang2019 fig. 3 but negative gradient in first layer
    TestNumberOfInterfaces() :
            model({{0_meter, 100_meter, 2600_meter_per_second},
                   {100_meter, 200_meter, 2400_meter_per_second},
                   {200_meter, 300_meter, 2400_meter_per_second},
                   {300_meter, 400_meter, 2700_meter_per_second},
                   {400_meter, 500_meter, 2250_meter_per_second}},
                  1000_meter, 1000_meter) {}

    VelocityModel model;
};

std::vector<DepthIntData> values = {{0_meter, 500_meter, 4},
                                    {50_meter, 450_meter, 4},
                                    {50_meter, 150_meter, 1},
                                    {200_meter, 250_meter, 0}};

INSTANTIATE_TEST_SUITE_P(TestInOrder, TestNumberOfInterfaces, testing::ValuesIn(values));


class TestModelGeometry : public testing::Test {
protected:
    TestModelGeometry() :
            model({{0_meter, 500_meter, 1_meter_per_second},
                   {5000_meter, 1000_meter, 1_meter_per_second}},
                  1000_meter, 1000_meter) {}

    VelocityModel model;
};

TEST_F(TestModelGeometry, TestAboveModel) {
    ASSERT_FALSE(model.in_model(0_meter, 0_meter, -2_meter));
}

TEST_F(TestModelGeometry, TestBelowModel) {
    ASSERT_FALSE(model.in_model(0_meter, 0_meter, 1001_meter));
}

class TestModelGeometryParametrized
        : public TestModelGeometry,
          public testing::WithParamInterface<std::tuple<Meter, Meter, Meter>> {};

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