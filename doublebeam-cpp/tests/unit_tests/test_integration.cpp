#include <gtest/gtest.h>

#include "raytracing_helpers.hpp"
#include "raytracing_types.hpp"

TEST(SnellsLaw, CompareWithPython) {
    Velocity v_before(2000);
    Velocity v_after(2200);
    Slowness p_in(0.7248451923249044_second_per_meter, 0.205219943364641_second_per_meter,
                  0.969149611692033_second_per_meter);
    Slowness p_out(0.7248451923249044_second_per_meter, 0.205219943364641_second_per_meter,
                   0.9691495893072383_second_per_meter);

    auto [px, py, pz] = snells_law(p_in, v_before, v_after, WaveType::Transmitted);
    EXPECT_FALSE(pz == p_in.pz);
    EXPECT_DOUBLE_EQ(px.get(), p_out.px.get());
    EXPECT_DOUBLE_EQ(py.get(), p_out.py.get());
    EXPECT_DOUBLE_EQ(pz.get(), p_out.pz.get());
}

TEST(SnellsLaw, ReflectedRay) {
    /**
     * For a horizontal interface and only one wave type (resulting in no velocity change
     * for the reflected wave at the interface in the same layer), only the sign of pz is
     * changed.
     */
    Slowness p_in(0.7248451923249044_second_per_meter, 0.205219943364641_second_per_meter,
                  0.969149611692033_second_per_meter);
    auto [px, py, pz] =
        snells_law(p_in, 1000_meter_per_second, 2000_meter_per_second, WaveType::Reflected);
    EXPECT_DOUBLE_EQ(px.get(), p_in.px.get());
    EXPECT_DOUBLE_EQ(py.get(), p_in.py.get());
    EXPECT_DOUBLE_EQ(pz.get(), -p_in.pz.get());
}