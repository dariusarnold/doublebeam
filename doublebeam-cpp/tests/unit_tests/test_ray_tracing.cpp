#include "model.hpp"
#include "printing.hpp"
#include "ray.hpp"
#include "raytracing.hpp"
#include "testing_utils.hpp"
#include "gtest/gtest.h"


class TestRayTracingBase : public ::testing::Test {
public:
    TestRayTracingBase() :
            vm{{{0_meter, 100_meter, 2000_meter_per_second},
                {100_meter, 200_meter, 2400_meter_per_second},
                {200_meter, 300_meter, 2650_meter_per_second},
                {300_meter, 400_meter, 2700_meter_per_second},
                {400_meter, 500_meter, 2925_meter_per_second}},
               3000_meter,
               3000_meter},
            krt(vm) {}

    VelocityModel vm;
    RayTracer krt;
};

TEST_F(TestRayTracingBase, TestStopAtCertainDepthConstantVelocityLayerDownwards) {
    // This ray ends in a linear velocity gradient layer and travels downwards
    auto initial_state = init_state(0_meter, 0_meter, 100_meter, vm, 0_rad, 0_rad);
    auto ray = krt.trace_ray(initial_state, "", 150_meter).value();
    auto [x, y, z] = ray.last_position();
    EXPECT_EQ(x.get(), 0);
    EXPECT_EQ(y.get(), 0);
    EXPECT_EQ(z.get(), 150);
}

TEST_F(TestRayTracingBase, TestStopAtCertainDepthConstantVelocityLayerUpwards) {
    InverseVelocity p(-1 / vm.eval_at(0_meter, 0_meter, 199_meter).value().get());
    // This ray ends in a linear velocity gradient layer
    auto initial_state =
        make_state(0_meter, 0_meter, 199_meter, InverseVelocity(0), InverseVelocity(0), p);
    auto ray = krt.trace_ray(initial_state, "", 150_meter).value();
    auto [x, y, z] = ray.last_position();
    // px is very small due to numerical inaccuracies calc
    EXPECT_DOUBLE_EQ(x.get(), 0);
    EXPECT_EQ(y, 0_meter);
    EXPECT_EQ(z, 150_meter);
}

TEST_F(TestRayTracingBase, TestStopAtCertainDepthLinearVelocityLayer) {
    // This ray ends in a linear velocity gradient layer
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, vm, 0_rad, 0_rad);
    Meter stop_depth(50);
    auto ray = krt.trace_ray(initial_state, "", stop_depth).value();
    auto [x, y, z] = ray.last_position();
    EXPECT_EQ(x, 0_meter);
    EXPECT_EQ(y, 0_meter);
    EXPECT_EQ(z, stop_depth);
}

TEST_F(TestRayTracingBase, TestStopAtCertainDepthOneLayerBottomToTop) {
    // This ray ends in a linear velocity gradient layer
    auto initial_state = init_state(1_meter, 1_meter, 99_meter, vm, radians(180_deg), 0_rad);
    auto ray = krt.trace_ray(initial_state, "", 50_meter).value();
    auto [x, y, z] = ray.last_position();
    EXPECT_TRUE(Close(x.get(), 1.));
    EXPECT_DOUBLE_EQ(y.get(), 1.);
    EXPECT_DOUBLE_EQ(z.get(), 50.);
}

TEST_F(TestRayTracingBase, TestStopAtCertainDepthMultiLayer) {
    // This ray ends in constant velocity layer
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, vm, 0_rad, 0_rad);
    auto ray = krt.trace_ray(initial_state, "T", 150_meter).value();
    auto [x, y, z] = ray.last_position();
    EXPECT_EQ(x.get(), 0);
    EXPECT_EQ(y.get(), 0);
    EXPECT_EQ(z.get(), 150);
}

TEST_F(TestRayTracingBase, TestStopAtCertainDepthMultiLayerWithNonVerticalRay) {
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, vm, radians(15_deg), 0_rad);
    auto ray = krt.trace_ray(initial_state, "T", 150_meter).value();
    auto [x, y, z] = ray.last_position();
    EXPECT_NE(x.get(), 0);
    EXPECT_EQ(y.get(), 0);
    EXPECT_EQ(z.get(), 150);
}

class TestRayTracing
        : public TestRayTracingBase,
          public ::testing::WithParamInterface<std::pair<std::array<double, 3>, double>> {};

TEST_P(TestRayTracing, TestCorrectEndpoint) {
    auto [slowness, endpoint] = GetParam();
    auto [px, py, pz] = slowness;
    RayState initial_state{Position{0_meter, 0_meter, 0_meter},
                           Slowness{InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)},
                           TravelTime{0_second}, Arclength{0_meter}};
    auto ray = krt.trace_ray(initial_state, "TTTT").value();
    auto last_traced_point = ray.last_position();
    EXPECT_TRUE(Close(last_traced_point.x.get(), endpoint, 1E-6));
}

TEST_P(TestRayTracing, TestArcLengthIncreasesContinuously) {
    // Test if arclength stays the same between an interface crossing and if it increases along
    // the ray.
    auto [slowness, _] = GetParam();
    auto [px, py, pz] = slowness;
    RayState initial_state{Position{0_meter, 0_meter, 0_meter},
                           Slowness{InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)},
                           TravelTime{0_second}, Arclength{0_meter}};
    auto ray = krt.trace_ray(initial_state, "TTTT").value();
    Arclength arclength = ray[0].begin().arclength;
    EXPECT_EQ(arclength.length.get(), 0) << "Initial arc length not zero.";
    for (const auto& segment : ray) {
        EXPECT_EQ(segment.begin().arclength, arclength)
            << "Arclength changed when crossing interface.";
        EXPECT_GE(segment.end().arclength.length.get(), arclength.length.get())
            << "Arclength did not increase in layer.";
        arclength = segment.end().arclength;
    }
}

using test_data_t = std::vector<std::pair<std::array<double, 3>, double>>;

// values pre calculated from two point ray tracing
test_data_t test_data = {{{0.00026509666908597633, 0, 0.00042393838707944381}, 469},
                         {{0.00032493968941590484, 0, 0.00038001868143855151}, 868},
                         {{0.00034109782958686665, 0, 0.00036558483372690507}, 2159},
                         {{0.00034130759537050534, 0, 0.00036538900550290698}, 2411}};

INSTANTIATE_TEST_SUITE_P(TestCorrectEndpoints, TestRayTracing, testing::ValuesIn(test_data));


TEST_F(TestRayTracing, LeaveTopOfModel) {
    RayState initial_state = make_state(0_meter, 0_meter, 0_meter, InverseVelocity(0),
                                        InverseVelocity(0), InverseVelocity(-0.005));
    EXPECT_THROW(krt.trace_ray(initial_state, "TRT"), std::runtime_error)
        << "Kinematic ray tracing not throwing when ray leaves top of model.";
    EXPECT_THROW(krt.trace_beam(initial_state, 1_meter, hertz_to_angular(1_hertz), "TRT"),
                 std::runtime_error)
        << "Dynamic ray tracing not throwing when beam leaves top of model";
}

TEST_F(TestRayTracing, LeaveBottomOfModel) {
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, vm, radians(10_deg), 0_rad);
    EXPECT_THROW(krt.trace_ray(initial_state, "TTTTTT"), std::runtime_error);
    EXPECT_THROW(krt.trace_beam(initial_state, 1_meter, hertz_to_angular(1_hertz), "TTTTTT"),
                 std::runtime_error);
}