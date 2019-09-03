#include "gtest/gtest.h"
#include "testing_utils.hpp"
#include "integration.hpp"
#include "model.h"
#include "printing.h"


class TestRayTracing : public ::testing::TestWithParam<std::pair<std::array<double, 3>, double>> {
protected:
    TestRayTracing() : vm{{
                                  {0, 100, 1800, 4},
                                  {100, 200, 2400, 0},
                                  {200, 300, 2400, 1},
                                  {300, 400, 2700, 0},
                                  {400, 500, 2250, 1.5}
                          }}, krt(vm) {}

    VelocityModel vm;
    KinematicRayTracer krt;
};

TEST_P(TestRayTracing, TestCorrectEndpoint) {
    auto[slowness, endpoint] = GetParam();
    auto [px, py, pz] = slowness;
    state_type initial_state{0, 0, 0, px, py, pz, 0};
    auto ray = krt.trace_ray(initial_state, "TTTTRTTTT", 1, 1);
    auto last_traced_point = ray.segments.back().data.back();
    EXPECT_TRUE(Close(last_traced_point[Index::X], endpoint));
}

using test_data_t = std::vector<std::pair<std::array<double, 3>, double>>;

test_data_t test_data = {{{0.000166674323178,  0., 0.0005299638150872}, 469},
                         {{0.0002545148149717, 0., 0.0004938260668176}, 868},
                         {{0.0003320005004714, 0., 0.0004454409534331}, 2159},
                         {{0.000333271179152,  0., 0.0004444910532905}, 2411}};

INSTANTIATE_TEST_SUITE_P(TestCorrectEndpoints, TestRayTracing, testing::ValuesIn(test_data));