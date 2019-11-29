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
    SourceReceiverSlowness{{469 / 2, 0, 500}, {469, 0, 0}, std::sin(math::radians(30)) / 3000},
    SourceReceiverSlowness{{868 / 2, 0, 500}, {868, 0, 0}, std::sin(math::radians(50)) / 3000},
    SourceReceiverSlowness{{2159 / 2, 0, 500}, {2159, 0, 0}, std::sin(math::radians(85)) / 3000}};

INSTANTIATE_TEST_SUITE_P(CompareWithFig9FromFang2019, TestTwoPointRayTracingWithData,
                         testing::ValuesIn(source_receiver_combinations));

std::vector<SourceReceiverSlowness> receiver_source_combinations = {
    SourceReceiverSlowness{{469, 0, 0}, {469 / 2, 0, 500}, -std::sin(math::radians(30)) / 3000},
    SourceReceiverSlowness{{868, 0, 0}, {868 / 2, 0, 500}, -std::sin(math::radians(50)) / 3000},
    SourceReceiverSlowness{{2159, 0, 0}, {2159 / 2, 0, 500}, -std::sin(math::radians(85)) / 3000}};
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
    auto [px, py, pz] = twopoint.trace({0, 0, 3000}, {2500, 2500, 0});
    EXPECT_NE(px, 0);
    EXPECT_NE(py, 0);
    EXPECT_DOUBLE_EQ(px, py);
    EXPECT_TRUE(pz < 0);
}

TEST_F(TwopointConstantVelocityModelWithOneLayer, TestConstantVelocityLayerWithStopDepth) {
    double source_x = 5000, source_y = 5000, source_z = 0;
    double target_x = 6200, target_y = 6200, target_z = 2350;
    auto [px, py, pz] =
        twopoint.trace({source_x, source_y, source_z}, {target_x, target_y, target_z});
    std::cerr << px << " " << py << " " << pz;
    std::cerr << math::degrees(math::angle(px, py, pz, 0, 0, 1)) << std::endl;
    EXPECT_DOUBLE_EQ(pz, 1.68897e-4);
    EXPECT_DOUBLE_EQ(px, 8.624537e-5);
    EXPECT_DOUBLE_EQ(py, 8.624537e-5);
    auto tracer = RayTracer(model);
    auto ray = tracer.trace_ray(make_state(source_x, source_y, source_z, px, py, pz), "", target_z);
    auto [end_point_x, end_point_y, end_point_z] = last_point(ray.value());
    EXPECT_TRUE(Close(end_point_x, target_x, 1e-5));
    EXPECT_TRUE(Close(end_point_y, target_y, 1e-5));
    EXPECT_DOUBLE_EQ(end_point_z, target_z);
}