#include <gtest/gtest.h>

#include "raytracing_helpers.hpp"
#include "raytracing_types.hpp"

TEST(SnellsLaw, CompareWithPython) {
    double_t v_before = 2000, v_after = 2200;
    double px_in = 0.7248451923249044, py_in = 0.205219943364641, pz_in = 0.969149611692033;
    double px_out = 0.7248451923249044, py_out = 0.205219943364641, pz_out = 0.9691495893072383;

    auto [px, py, pz] = snells_law(px_in, py_in, pz_in, v_before, v_after, WaveType::Transmitted);
    EXPECT_FALSE(pz == pz_in);
    EXPECT_DOUBLE_EQ(px, px_out);
    EXPECT_DOUBLE_EQ(py, py_out);
    EXPECT_DOUBLE_EQ(pz, pz_out);
}

TEST(SnellsLaw, ReflectedRay) {
    /**
     * For a horizontal interface and only one wave type (resulting in no velocity change
     * for the reflected wave at the interface in the same layer), only the sign of pz is
     * changed.
     */
    double px_in = 0.7248451923249044, py_in = 0.205219943364641, pz_in = 0.969149611692033;
    auto [px, py, pz] = snells_law(px_in, py_in, pz_in, 1000, 2000, WaveType::Reflected);
    EXPECT_DOUBLE_EQ(px, px_in);
    EXPECT_DOUBLE_EQ(py, py_in);
    EXPECT_DOUBLE_EQ(pz, -pz_in);
}