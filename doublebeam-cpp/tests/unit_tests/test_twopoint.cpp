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
    EXPECT_THROW(twopoint.trace({0, 0, top - 1}, {0, 0, bottom}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, bottom + 1}, {0, 0, bottom}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, top}, {0, 0, top - 1}), std::domain_error);
    EXPECT_THROW(twopoint.trace({0, 0, top}, {0, 0, bottom + 1}), std::domain_error);
}


/**
 * Simple struct describing one source/receiver pair and the expected
 * slowness for two point ray tracing.
 */
struct SourceReceiverSlowness {
    position_t source;
    position_t receiver;
    double px;
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
    EXPECT_TRUE(Close(px, px_target, 1e-5, 2e-6))
        << "For source " << source << ", receiver " << receiver;
}

std::vector<SourceReceiverSlowness> source_receiver_combinations = {
    SourceReceiverSlowness{{469. / 2, 0, 500}, {469, 0, 0}, std::sin(math::radians(30)) / 3000},
    SourceReceiverSlowness{{868. / 2, 0, 500}, {868, 0, 0}, std::sin(math::radians(50)) / 3000},
    SourceReceiverSlowness{{2159. / 2, 0, 500}, {2159, 0, 0}, std::sin(math::radians(85)) / 3000}};

INSTANTIATE_TEST_SUITE_P(CompareWithFig9FromFang2019, TestTwoPointRayTracingWithData,
                         testing::ValuesIn(source_receiver_combinations));

std::vector<SourceReceiverSlowness> receiver_source_combinations = {
    SourceReceiverSlowness{{469, 0, 0}, {469. / 2, 0, 500}, -std::sin(math::radians(30)) / 3000},
    SourceReceiverSlowness{{868, 0, 0}, {868. / 2, 0, 500}, -std::sin(math::radians(50)) / 3000},
    SourceReceiverSlowness{{2159, 0, 0}, {2159. / 2, 0, 500}, -std::sin(math::radians(85)) / 3000}};
INSTANTIATE_TEST_SUITE_P(SwapSourceAndReceiver, TestTwoPointRayTracingWithData,
                         testing::ValuesIn(receiver_source_combinations));

INSTANTIATE_TEST_SUITE_P(
    3DCase, TestTwoPointRayTracingWithData,
    testing::Values(SourceReceiverSlowness{{0, 0, 500}, {100, 200, 0}, 7.16965157e-05},
                    SourceReceiverSlowness{{100, 200, 0}, {0, 0, 500}, -7.16965157e-05}));

TEST_F(TestTwopointRayTracingBase, TestStraightDown) {
    auto [px, py, pz] = twopoint.trace({0, 0, 0}, {0, 0, 400});
    EXPECT_EQ(px, 0);
    EXPECT_EQ(py, 0);
    EXPECT_DOUBLE_EQ(pz, 1 / model.eval_at(0, 0, 0).value());
}

TEST_F(TestTwopointRayTracingBase, TestStraightDownWithXYCoordinates) {
    auto [px, py, pz] = twopoint.trace({100, 100, 450}, {100, 100, 0});
    EXPECT_EQ(px, 0);
    EXPECT_EQ(py, 0);
    EXPECT_EQ(pz, -1. / model.eval_at(100, 100, 450).value());
}

TEST_F(TestTwopointRayTracingBase, TestNanReturn) {
    // This test fails due to a Nan which has to be summed using Nansum. Check if code handles
    // this case correctly.
    auto [px, py, pz] = twopoint.trace({400, 500, 0}, {400, 400, 450});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
}

// This works
TEST_F(TestTwopointRayTracingBase, Test3DCaseBottomToTop) {
    auto [px, py, pz] = twopoint.trace({25, 50, 450}, {100, 100, 0});
    auto px_expected = 6.56403315e-05, py_expected = 4.37602210e-05, pz_expected = -3.32653811e-04;
    EXPECT_TRUE(Close(px, px_expected, 1e-5, 2e-6));
    EXPECT_TRUE(Close(py, py_expected, 1e-5, 2e-6));
    EXPECT_TRUE(Close(pz, pz_expected, 1e-5, 2e-6));
}

// This works with the workaround but not without
TEST_F(TestTwopointRayTracingBase, Test3DCaseTopToBottom) {
    auto [px, py, pz] = twopoint.trace({100, 100, 0}, {25, 50, 450});
    auto px_expected = -5.81890528e-05, py_expected = -3.87927019e-05, pz_expected = 5.51136222e-04;
    EXPECT_TRUE(Close(px, px_expected, 1e-5, 6e-5));
    EXPECT_TRUE(Close(py, py_expected, 1e-5, 6e-5));
    EXPECT_TRUE(Close(pz, pz_expected, 1e-5, 6e-5));
}


TEST_F(TestTwopointRayTracingBase, TestNaNIsFixedTopToBottom) {
    // This used to give an error due to omega_k_tilde being negative in a sqrt since omega_k
    // (the velocity at the bottom of a layer) was larger than v_M (the largest velocity
    auto [px, py, pz] = twopoint.trace({5000, 5000, 0}, {6200, 6200, 450});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
    EXPECT_DOUBLE_EQ(px, py);
    EXPECT_NE(px, 0);
}

TEST_F(TestTwopointRayTracingBase, TestNaNIsFixedBottomToTop) {
    // This used to give an error due to omega_k_tilde being negative in a sqrt since omega_k
    // (the velocity at the bottom of a layer) was larger than v_M (the largest velocity
    auto [px, py, pz] = twopoint.trace({6200, 6200, 450}, {5000, 5000, 0});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
    EXPECT_DOUBLE_EQ(px, py);
    EXPECT_NE(px, 0);
}

TEST_F(TestTwopointRayTracingBase, TestBottomToTopEquivalencyWith2D) {
    auto [px, py, pz] = twopoint.trace({0, 0, 450}, {1200, 0, 0});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
    // this result was obtained using the Matlab code
    EXPECT_TRUE(Close(px, 3.3284e-04, 1e-7, 1e-5));
    EXPECT_DOUBLE_EQ(py, 0);
}

TEST_F(TestTwopointRayTracingBase, TestInSingleLayerBottomToTop) {
    auto [px, py, pz] = twopoint.trace({125, 0, 50}, {0, 0, 0});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
    EXPECT_NE(px, 0);
    EXPECT_TRUE(Close(py, 0., 0., 1e-19));
}

TEST_F(TestTwopointRayTracingBase, TestInSingleLayerTopToBottom) {
    auto [px, py, pz] = twopoint.trace({0, 0, 0}, {125, 0, 50});
    EXPECT_FALSE(std::isnan(px));
    EXPECT_FALSE(std::isnan(py));
    EXPECT_FALSE(std::isnan(pz));
    EXPECT_NE(px, 0);
    EXPECT_DOUBLE_EQ(py, 0);
    EXPECT_NE(pz, 0);
}

class TwopointConstantVelocityModelWithOneLayer : public testing::Test {
protected:
    double constant_velocity = 4800;
    VelocityModel model{{Layer{0, 3000, constant_velocity, 0}}, 10000, 10000};
    TwoPointRayTracing twopoint{model};
};

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestVerticalRayDownwards) {
    auto [px, py, pz] = twopoint.trace({0, 0, 0}, {0, 0, 2500});
    double expected_pz = 1 / constant_velocity;
    EXPECT_DOUBLE_EQ(px, 0);
    EXPECT_DOUBLE_EQ(py, 0);
    EXPECT_DOUBLE_EQ(pz, expected_pz);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestVerticalRayUpwards) {
    auto [px, py, pz] = twopoint.trace({0, 0, 2500}, {0, 0, 0});
    double expected_pz = -1 / constant_velocity;
    EXPECT_DOUBLE_EQ(px, 0);
    EXPECT_DOUBLE_EQ(py, 0);
    EXPECT_DOUBLE_EQ(pz, expected_pz);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestLongDistance) {
    double sx = 0, sy = 0, sz = 3000;
    double rx = 2500, ry = 2500, rz = 0;
    auto [px, py, pz] = twopoint.trace({sx, sy, sz}, {rx, ry, rz});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(sx, sy, sz, px, py, pz));
    auto [x, y, z] = last_point(ray.value());
    EXPECT_TRUE(Close(x, rx));
    EXPECT_TRUE(Close(y, ry));
    EXPECT_TRUE(Close(z, rz));
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestLongDistanceSwapped) {
    double sx = 2500, sy = 2500, sz = 0;
    double rx = 0, ry = 0, rz = 3000;
    auto [px, py, pz] = twopoint.trace({sx, sy, sz}, {rx, ry, rz});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(sx, sy, sz, px, py, pz));
    auto [x, y, z] = last_point(ray.value());
    EXPECT_TRUE(Close(x, rx, 0., 1e-12));
    EXPECT_TRUE(Close(y, ry, 0., 1e-12));
    EXPECT_TRUE(Close(z, rz));
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestConstantVelocityLayerWithStopDepth) {
    double source_x = 5000, source_y = 5000, source_z = 0;
    double target_x = 6200, target_y = 6200, target_z = 2350;
    auto [px, py, pz] =
        twopoint.trace({source_x, source_y, source_z}, {target_x, target_y, target_z});
    std::cerr << px << " " << py << " " << pz;
    std::cerr << math::degrees(math::angle(px, py, pz, 0, 0, 1)) << std::endl;
    EXPECT_TRUE(Close(px, 8.624537e-5));
    EXPECT_TRUE(Close(py, 8.624537e-5));
    EXPECT_TRUE(Close(pz, 1.68897e-4, 1e-6));
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(source_x, source_y, source_z, px, py, pz), "", target_z);
    auto [end_point_x, end_point_y, end_point_z] = last_point(ray.value());
    EXPECT_TRUE(Close(end_point_x, target_x, 1e-5));
    EXPECT_TRUE(Close(end_point_y, target_y, 1e-5));
    EXPECT_DOUBLE_EQ(end_point_z, target_z);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestConstantVelocityLayerWithStopDepthUpwards) {
    double source_x = 6200, source_y = 6200, source_z = 2350;
    double target_x = 5000, target_y = 5000, target_z = 0;
    auto [px, py, pz] =
        twopoint.trace({source_x, source_y, source_z}, {target_x, target_y, target_z});
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(source_x, source_y, source_z, px, py, pz), "", target_z);
    auto [end_point_x, end_point_y, end_point_z] = last_point(ray.value());
    EXPECT_TRUE(Close(end_point_x, target_x, 1e-5));
    EXPECT_TRUE(Close(end_point_y, target_y, 1e-5));
    EXPECT_DOUBLE_EQ(end_point_z, target_z);
}


/**
 * Test two point ray tracing with 10 pregenerated random points.
 * The slowness is used to trace a ray and then the endpoint of the ray is compared with the target.
 * This case uses a constant velocity model with a single layer.
 */
class TestTwoPointRayTracingOneLayer : public TwopointConstantVelocityModelWithOneLayer,
                                       public testing::WithParamInterface<position_t> {
protected:
    RayTracer ray_tracer{model};
};


TEST_P(TestTwoPointRayTracingOneLayer, TestRandomCases) {
    auto start = position_t{5000, 5000, 1500};
    position_t target = GetParam();
    auto slowness = twopoint.trace(start, target);
    auto ray = ray_tracer.trace_ray(make_state(start, slowness), "", std::get<2>(target));
    auto [x, y, z] = last_point(ray.value());
    auto [target_x, target_y, target_z] = target;
    EXPECT_TRUE(Close(x, target_x, 1e-5));
    EXPECT_TRUE(Close(y, target_y, 1e-5));
    EXPECT_DOUBLE_EQ(z, target_z);
}

std::vector<position_t> get_points() {
    // return 10 random points between 0 and 10000
    std::vector<double> x_coords{2999.93680517, 5673.46859475, 2901.02666565, 2239.096261,
                                 6631.67512127, 8086.75945027, 2563.18131644, 7944.76077953,
                                 259.91226483,  1292.27259846};

    std::vector<double> y_coords{590.56991678,  9136.3913014, 3103.69739044, 2282.05573332,
                                 5538.70658045, 813.79285535, 3554.28191681, 1206.13755077,
                                 5431.61049355, 8203.57325616};

    std::vector<double> z_coords{70.77988368,   2373.77000553, 583.790722,   2214.03316274,
                                 1127.69456278, 2687.38474072, 915.10591529, 1390.67106024,
                                 2793.65525894, 2068.9689757};

    std::vector<position_t> points;
    int i = -1;
    std::generate_n(std::back_inserter(points), z_coords.size(), [&]() {
        i++;
        return position_t{x_coords[i], y_coords[i], z_coords[i]};
    });
    return points;
}

const std::vector<position_t> points_single_layer = get_points();

INSTANTIATE_TEST_SUITE_P(TestManyRandomPoints, TestTwoPointRayTracingOneLayer,
                         testing::ValuesIn(points_single_layer));


class TestTwoPointRayTracingMultipleLayer : public TestTwopointRayTracingBase,
                                            public testing::WithParamInterface<position_t> {
protected:
    RayTracer ray_tracer{model};
};

std::vector<position_t> scale_z_coordinate() {
    std::vector<position_t> old_points = get_points();
    double new_max_depth = 500;
    double old_max_depth = 10000;
    for (auto& p : old_points) {
        std::get<2>(p) *= new_max_depth / old_max_depth;
    }
    return old_points;
}

std::vector<WaveType> direct_ray_code(position_t source, position_t receiver,
                                      const VelocityModel& model) {
    auto n = model.number_of_interfaces_between(std::get<Index::Z>(source),
                                                std::get<Index::Z>(receiver));
    return std::vector<WaveType>(n, WaveType::Transmitted);
}

TEST_P(TestTwoPointRayTracingMultipleLayer, TestByComparisionWithRayTracing) {
    auto start = position_t{5000, 5000, 250};
    position_t target = GetParam();
    auto slowness = twopoint.trace(start, target);
    auto ray_code = direct_ray_code(start, target, model);
    auto ray = ray_tracer.trace_ray(make_state(start, slowness), ray_code, std::get<2>(target));
    auto [x, y, z] = last_point(ray.value());
    auto [target_x, target_y, target_z] = target;
    EXPECT_TRUE(Close(x, target_x, 1e-5));
    EXPECT_TRUE(Close(y, target_y, 1e-5));
    EXPECT_DOUBLE_EQ(z, target_z);
}

const std::vector<position_t> points_multi_layer = scale_z_coordinate();

INSTANTIATE_TEST_SUITE_P(TestManyRandomPoints, TestTwoPointRayTracingMultipleLayer,
                         testing::ValuesIn(points_multi_layer));