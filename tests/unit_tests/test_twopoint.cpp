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
#include "gtest/gtest.h"
#include <testing_utils.hpp>

#include "model.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iterator>


class TestTwopointRayTracingBase : public ::testing::Test {
protected:
    TestTwopointRayTracingBase() : twopoint(model) {}
    VelocityModel model =
        read_velocity_file(current_source_path(__FILE__) / "data" / "twopoint_velocity_model.txt");
    TwoPointRayTracing twopoint;
};


/**
 * Test if sources/receivers outside of the velocity model raise an exception.
 */
TEST_F(TestTwopointRayTracingBase, TestThrowWhenOutOfRange) {
    auto [top, bottom] = model.get_top_bottom();
    EXPECT_THROW(twopoint.trace({0_meter, 0_meter, top - 1_meter}, {0_meter, 0_meter, bottom}),
                 std::domain_error);
    EXPECT_THROW(twopoint.trace({0_meter, 0_meter, bottom + 1_meter}, {0_meter, 0_meter, bottom}),
                 std::domain_error);
    EXPECT_THROW(twopoint.trace({0_meter, 0_meter, top}, {0_meter, 0_meter, top - 1_meter}),
                 std::domain_error);
    EXPECT_THROW(twopoint.trace({0_meter, 0_meter, top}, {0_meter, 0_meter, bottom + 1_meter}),
                 std::domain_error);
}


/**
 * Simple struct describing one source/receiver pair and the expected
 * slowness for two point ray tracing.
 */
struct SourceReceiverSlowness {
    Position source;
    Position receiver;
    InverseVelocity px;
};

std::ostream& operator<<(std::ostream& os, std::tuple<double, double, double> t) {
    os << "(" << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, SourceReceiverSlowness d) {
    os << "Source: " << d.source << " Receiver: " << d.receiver << " Slowness: " << d.px;
    return os;
}

class TestTwoPointRayTracingWithData
        : public TestTwopointRayTracingBase,
          public ::testing::WithParamInterface<SourceReceiverSlowness> {};

/**
 * Test if two point ray tracing gives the same slowness as calculated from fig. 9 in Fang2019.
 * The velocity model has all linear velocity layers replaced by a constant velocity layer with the
 * average velocity. This way results are only slightly inaccurate and can be used as a good first
 * order test if the method is working.
 */
TEST_P(TestTwoPointRayTracingWithData, CompareWithGivenResult) {
    auto [source, receiver, px_target] = GetParam();
    auto [px, py, pz] = twopoint.trace(source, receiver);
    EXPECT_TRUE(Close(px.get(), px_target.get(), 1e-5, 2e-6))
        << "For source " << source << ", receiver " << receiver;
}

std::vector<SourceReceiverSlowness> source_receiver_combinations = {
    SourceReceiverSlowness{{469_meter / 2, 0_meter, 500_meter},
                           {469_meter, 0_meter, 0_meter},
                           InverseVelocity(std::sin(radians(30_deg).get()) / 3000)},
    SourceReceiverSlowness{{868_meter / 2, 0_meter, 500_meter},
                           {868_meter, 0_meter, 0_meter},
                           InverseVelocity(std::sin(radians(50_deg).get()) / 3000)},
    SourceReceiverSlowness{{2159_meter / 2, 0_meter, 500_meter},
                           {2159_meter, 0_meter, 0_meter},
                           InverseVelocity(std::sin(radians(85_deg).get()) / 3000)}};

INSTANTIATE_TEST_SUITE_P(CompareWithFig9FromFang2019, TestTwoPointRayTracingWithData,
                         testing::ValuesIn(source_receiver_combinations));

// swap source and receiver position
std::vector<SourceReceiverSlowness>
    receiver_source_combinations(source_receiver_combinations.rbegin(),
                                 source_receiver_combinations.rend());
INSTANTIATE_TEST_SUITE_P(SwapSourceAndReceiver, TestTwoPointRayTracingWithData,
                         testing::ValuesIn(receiver_source_combinations));

INSTANTIATE_TEST_SUITE_P(3DCase, TestTwoPointRayTracingWithData,
                         testing::Values(SourceReceiverSlowness{{0_meter, 0_meter, 500_meter},
                                                                {100_meter, 200_meter, 0_meter},
                                                                7.16965157e-05_second_per_meter},
                                         SourceReceiverSlowness{{100_meter, 200_meter, 0_meter},
                                                                {0_meter, 0_meter, 500_meter},
                                                                -7.16965157e-05_second_per_meter}));

TEST_F(TestTwopointRayTracingBase, TestStraightDown) {
    auto [px, py, pz] = twopoint.trace({0_meter, 0_meter, 0_meter}, {0_meter, 0_meter, 400_meter});
    EXPECT_EQ(px, 0_second_per_meter);
    EXPECT_EQ(py, 0_second_per_meter);
    EXPECT_DOUBLE_EQ(pz.get(), 1 / model.eval_at(0_meter, 0_meter, 0_meter).value().get());
}

TEST_F(TestTwopointRayTracingBase, TestStraightDownWithXYCoordinates) {
    auto [px, py, pz] =
        twopoint.trace({100_meter, 100_meter, 450_meter}, {100_meter, 100_meter, 0_meter});
    EXPECT_EQ(px, 0_second_per_meter);
    EXPECT_EQ(py, 0_second_per_meter);
    EXPECT_EQ(pz.get(), -1. / model.eval_at(100_meter, 100_meter, 450_meter).value().get());
}

TEST_F(TestTwopointRayTracingBase, TestNanReturn) {
    // This test fails due to a Nan which has to be summed using Nansum. Check if code handles
    // this case correctly.
    auto [px, py, pz] =
        twopoint.trace({400_meter, 500_meter, 0_meter}, {400_meter, 400_meter, 450_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
}

// This works
TEST_F(TestTwopointRayTracingBase, Test3DCaseBottomToTop) {
    auto [px, py, pz] =
        twopoint.trace({25_meter, 50_meter, 450_meter}, {100_meter, 100_meter, 0_meter});
    auto px_expected = 6.56403315e-05, py_expected = 4.37602210e-05, pz_expected = -3.32653811e-04;
    EXPECT_TRUE(Close(px.get(), px_expected, 1e-5, 2e-6));
    EXPECT_TRUE(Close(py.get(), py_expected, 1e-5, 2e-6));
    EXPECT_TRUE(Close(pz.get(), pz_expected, 1e-5, 2e-6));
}

// This works with the workaround but not without
TEST_F(TestTwopointRayTracingBase, Test3DCaseTopToBottom) {
    auto [px, py, pz] =
        twopoint.trace({100_meter, 100_meter, 0_meter}, {25_meter, 50_meter, 450_meter});
    auto px_expected = -5.81890528e-05, py_expected = -3.87927019e-05, pz_expected = 5.51136222e-04;
    EXPECT_TRUE(Close(px.get(), px_expected, 1e-5, 6e-5));
    EXPECT_TRUE(Close(py.get(), py_expected, 1e-5, 6e-5));
    EXPECT_TRUE(Close(pz.get(), pz_expected, 1e-5, 6e-5));
}


TEST_F(TestTwopointRayTracingBase, TestNaNIsFixedTopToBottom) {
    // This used to give an error due to omega_k_tilde being negative in a sqrt since omega_k
    // (the velocity at the bottom of a layer) was larger than v_M (the largest velocity
    auto [px, py, pz] =
        twopoint.trace({5000_meter, 5000_meter, 0_meter}, {6200_meter, 6200_meter, 450_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
    EXPECT_DOUBLE_EQ(px.get(), py.get());
    EXPECT_NE(px, 0_second_per_meter);
}

TEST_F(TestTwopointRayTracingBase, TestNaNIsFixedBottomToTop) {
    // This used to give an error due to omega_k_tilde being negative in a sqrt since omega_k
    // (the velocity at the bottom of a layer) was larger than v_M (the largest velocity
    auto [px, py, pz] =
        twopoint.trace({6200_meter, 6200_meter, 450_meter}, {5000_meter, 5000_meter, 0_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
    EXPECT_DOUBLE_EQ(px.get(), py.get());
    EXPECT_NE(px, 0_second_per_meter);
}

TEST_F(TestTwopointRayTracingBase, TestBottomToTopEquivalencyWith2D) {
    auto [px, py, pz] =
        twopoint.trace({0_meter, 0_meter, 450_meter}, {1200_meter, 0_meter, 0_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
    // this result was obtained using the Matlab code
    EXPECT_TRUE(Close(px.get(), 3.3284e-04, 1e-7, 1e-5));
    EXPECT_DOUBLE_EQ(py.get(), 0);
}

TEST_F(TestTwopointRayTracingBase, TestInSingleLayerBottomToTop) {
    auto [px, py, pz] = twopoint.trace({125_meter, 0_meter, 50_meter}, {0_meter, 0_meter, 0_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
    EXPECT_NE(px, 0_second_per_meter);
    EXPECT_TRUE(Close(py.get(), 0., 0., 1e-19));
}

TEST_F(TestTwopointRayTracingBase, TestInSingleLayerTopToBottom) {
    auto [px, py, pz] = twopoint.trace({0_meter, 0_meter, 0_meter}, {125_meter, 0_meter, 50_meter});
    EXPECT_FALSE(isnan(px));
    EXPECT_FALSE(isnan(py));
    EXPECT_FALSE(isnan(pz));
    EXPECT_NE(px, 0_second_per_meter);
    EXPECT_DOUBLE_EQ(py.get(), 0);
    EXPECT_NE(pz, 0_second_per_meter);
}

class TwopointConstantVelocityModelWithOneLayer : public testing::Test {
protected:
    Velocity constant_velocity{4800};
    VelocityModel model{{Layer{0_meter, 3000_meter, constant_velocity}}, 10000_meter, 10000_meter};
    TwoPointRayTracing twopoint{model};
};

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestVerticalRayDownwards) {
    auto [px, py, pz] = twopoint.trace({0_meter, 0_meter, 0_meter}, {0_meter, 0_meter, 2500_meter});
    double expected_pz = 1 / constant_velocity.get();
    EXPECT_DOUBLE_EQ(px.get(), 0);
    EXPECT_DOUBLE_EQ(py.get(), 0);
    EXPECT_DOUBLE_EQ(pz.get(), expected_pz);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestVerticalRayUpwards) {
    auto [px, py, pz] = twopoint.trace({0_meter, 0_meter, 2500_meter}, {0_meter, 0_meter, 0_meter});
    double expected_pz = -1 / constant_velocity.get();
    EXPECT_DOUBLE_EQ(px.get(), 0);
    EXPECT_DOUBLE_EQ(py.get(), 0);
    EXPECT_DOUBLE_EQ(pz.get(), expected_pz);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestLongDistance) {
    auto sx = 0_meter, sy = 0_meter, sz = 3000_meter;
    auto rx = 2500_meter, ry = 2500_meter, rz = 0_meter;
    auto [px, py, pz] = twopoint.trace({sx, sy, sz}, {rx, ry, rz});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(
        make_state(sx, sy, sz, InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)));
    auto [x, y, z] = ray.value().last_position();
    EXPECT_TRUE(Close(x.get(), rx.get()));
    EXPECT_TRUE(Close(y.get(), ry.get()));
    EXPECT_TRUE(Close(z.get(), rz.get()));
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestLongDistanceSwapped) {
    auto sx = 2500_meter, sy = 2500_meter, sz = 0_meter;
    auto rx = 0_meter, ry = 0_meter, rz = 3000_meter;
    auto [px, py, pz] = twopoint.trace({sx, sy, sz}, {rx, ry, rz});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(
        make_state(sx, sy, sz, InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)));
    auto [x, y, z] = ray.value().last_position();
    EXPECT_TRUE(Close(x.get(), rx.get(), 0., 1e-12));
    EXPECT_TRUE(Close(y.get(), ry.get(), 0., 1e-12));
    EXPECT_TRUE(Close(z.get(), rz.get()));
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestConstantVelocityLayerWithStopDepth) {
    auto source_x = 5000_meter, source_y = 5000_meter, source_z = 0_meter;
    auto target_x = 6200_meter, target_y = 6200_meter, target_z = 2350_meter;
    auto [px, py, pz] =
        twopoint.trace({source_x, source_y, source_z}, {target_x, target_y, target_z});
    EXPECT_TRUE(Close(px.get(), 8.624537e-5));
    EXPECT_TRUE(Close(py.get(), 8.624537e-5));
    EXPECT_TRUE(Close(pz.get(), 1.68897e-4, 1e-6));
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(source_x, source_y, source_z, InverseVelocity(px),
                                           InverseVelocity(py), InverseVelocity(pz)),
                                "", target_z);
    auto [end_point_x, end_point_y, end_point_z] = ray.value().last_position();
    EXPECT_TRUE(Close(end_point_x.get(), target_x.get(), 1e-5));
    EXPECT_TRUE(Close(end_point_y.get(), target_y.get(), 1e-5));
    EXPECT_DOUBLE_EQ(end_point_z.get(), target_z.get());
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestConstantVelocityLayerWithStopDepthUpwards) {
    auto source_x = 6200_meter, source_y = 6200_meter, source_z = 2350_meter;
    auto target_x = 5000_meter, target_y = 5000_meter, target_z = 0_meter;
    auto [px, py, pz] =
        twopoint.trace({source_x, source_y, source_z}, {target_x, target_y, target_z});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(source_x, source_y, source_z, InverseVelocity(px),
                                           InverseVelocity(py), InverseVelocity(pz)),
                                "", target_z);
    auto [end_point_x, end_point_y, end_point_z] = ray.value().last_position();
    EXPECT_TRUE(Close(end_point_x.get(), target_x.get(), 1e-5));
    EXPECT_TRUE(Close(end_point_y.get(), target_y.get(), 1e-5));
    EXPECT_DOUBLE_EQ(end_point_z.get(), target_z.get());
}


/**
 * Test two point ray tracing with 10 pregenerated random points.
 * The slowness is used to trace a ray and then the endpoint of the ray is compared with the target.
 * This case uses a constant velocity model with a single layer.
 */
class TestTwoPointRayTracingOneLayer : public TwopointConstantVelocityModelWithOneLayer,
                                       public testing::WithParamInterface<Position> {
protected:
    RayTracer ray_tracer{model};
};


TEST_P(TestTwoPointRayTracingOneLayer, TestRandomCases) {
    Position start{5000_meter, 5000_meter, 1500_meter};
    // TODO change twopoint api to accept Position instead of position_t
    Position start_pos(5000_meter, 5000_meter, 1500_meter);
    Position target = GetParam();
    auto [px, py, pz] = twopoint.trace(start, target);
    auto ray = ray_tracer.trace_ray(
        make_state(start_pos,
                   Slowness{InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)}),
        "", target.z);
    auto [x, y, z] = ray.value().last_position();
    auto [target_x, target_y, target_z] = target;
    EXPECT_TRUE(Close(x.get(), target_x.get(), 1e-5));
    EXPECT_TRUE(Close(y.get(), target_y.get(), 1e-5));
    EXPECT_DOUBLE_EQ(z.get(), target_z.get());
}

std::vector<Position> get_points() {
    // return 10 random points between 0 and 10000
    std::vector<Meter> x_coords{2999.93680517_meter, 5673.46859475_meter, 2901.02666565_meter,
                                2239.096261_meter,   6631.67512127_meter, 8086.75945027_meter,
                                2563.18131644_meter, 7944.76077953_meter, 259.91226483_meter,
                                1292.27259846_meter};

    std::vector<Meter> y_coords{590.56991678_meter,  9136.3913014_meter,  3103.69739044_meter,
                                2282.05573332_meter, 5538.70658045_meter, 813.79285535_meter,
                                3554.28191681_meter, 1206.13755077_meter, 5431.61049355_meter,
                                8203.57325616_meter};

    std::vector<Meter> z_coords{70.77988368_meter,   2373.77000553_meter, 583.790722_meter,
                                2214.03316274_meter, 1127.69456278_meter, 2687.38474072_meter,
                                915.10591529_meter,  1390.67106024_meter, 2793.65525894_meter,
                                2068.9689757_meter};

    std::vector<Position> points;
    int i = -1;
    std::generate_n(std::back_inserter(points), z_coords.size(), [&]() {
        i++;
        return Position{x_coords[i], y_coords[i], z_coords[i]};
    });
    return points;
}

const std::vector<Position> points_single_layer = get_points();

INSTANTIATE_TEST_SUITE_P(TestManyRandomPoints, TestTwoPointRayTracingOneLayer,
                         testing::ValuesIn(points_single_layer));


class TestTwoPointRayTracingMultipleLayer : public TestTwopointRayTracingBase,
                                            public testing::WithParamInterface<Position> {
protected:
    RayTracer ray_tracer{model};
};

std::vector<Position> scale_z_coordinate() {
    std::vector<Position> old_points = get_points();
    double new_max_depth = 500;
    double old_max_depth = 10000;
    for (auto& p : old_points) {
        p.z *= new_max_depth / old_max_depth;
    }
    return old_points;
}

std::vector<WaveType> direct_ray_code(Position source, Position receiver,
                                      const VelocityModel& model) {
    auto n = model.number_of_interfaces_between(source.z, receiver.z);
    return std::vector<WaveType>(n, WaveType::Transmitted);
}

TEST_P(TestTwoPointRayTracingMultipleLayer, TestByComparisionWithRayTracing) {
    auto start = Position{5000_meter, 5000_meter, 250_meter};
    Position start_pos(5000_meter, 5000_meter, 250_meter);
    Position target = GetParam();
    auto [px, py, pz] = twopoint.trace(start, target);
    auto ray_code = direct_ray_code(start, target, model);
    auto ray = ray_tracer.trace_ray(
        make_state(start_pos,
                   Slowness(InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz))),
        ray_code, target.z);
    auto [x, y, z] = ray.value().last_position();
    auto [target_x, target_y, target_z] = target;
    EXPECT_TRUE(Close(x.get(), target_x.get(), 1e-5));
    EXPECT_TRUE(Close(y.get(), target_y.get(), 1e-5));
    EXPECT_DOUBLE_EQ(z.get(), target_z.get());
}

const std::vector<Position> points_multi_layer = scale_z_coordinate();

INSTANTIATE_TEST_SUITE_P(TestManyRandomPoints, TestTwoPointRayTracingMultipleLayer,
                         testing::ValuesIn(points_multi_layer));