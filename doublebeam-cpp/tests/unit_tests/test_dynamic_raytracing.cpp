#include <gtest/gtest.h>

#include "model.hpp"
#include "dynamic_raytracing.hpp"


class DynamicRaytracingBase : public testing::Test {
protected:
    VelocityModel model = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    DynamicRayTracer rt{model};
};

TEST_F(DynamicRaytracingBase, TestThrowWhenStartingOutOfModel) {
    auto [top, bottom] = model.get_top_bottom();
    auto above = init_state(0, 0, top - 1, model, 0, 0, 0);
    auto below = init_state(0, 0, bottom + 1, model, 0, 0, 0);
    EXPECT_THROW(rt.trace_beam(above, 0, 0), std::domain_error);
    EXPECT_THROW(rt.trace_beam(below, 0, 0), std::domain_error);
}
