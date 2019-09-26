#include <gtest/gtest.h>

#include <xtensor/xio.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xview.hpp>

#include "kinematic_raytracing.hpp"
#include "model.hpp"
#include "utils.hpp"


class DynamicRaytracingBase : public testing::Test {
protected:
    VelocityModel model = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    KinematicRayTracer rt{model};
};

TEST_F(DynamicRaytracingBase, TestThrowWhenStartingOutOfModel) {
    auto [top, bottom] = model.get_top_bottom();
    auto above = init_state(0, 0, top - 1, model, 0, 0, 0);
    auto below = init_state(0, 0, bottom + 1, model, 0, 0, 0);
    EXPECT_THROW(rt.trace_beam(above, 0, 0), std::domain_error);
    EXPECT_THROW(rt.trace_beam(below, 0, 0), std::domain_error);
}

TEST_F(DynamicRaytracingBase, TestBasicCredibilityOfResult) {
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto width = 1;
    auto frequency = 2;
    auto beam = rt.trace_beam(initial_state, width, frequency, "TTTT");
    EXPECT_EQ(beam.segments.size(), 5) << "Not matching expected number of segments.";
    EXPECT_EQ(beam.width(), width) << "Beam width not set correctly.";
    EXPECT_EQ(beam.frequency(), frequency) << "Beam frequency not set correctly.";
    for (auto i = 0; i < beam.segments.size(); ++i) {
        EXPECT_EQ(beam[i].P.shape()[0], beam[i].Q.shape()[0]) << "Length of P and Q different.";
        EXPECT_EQ(beam[i].P.shape()[0], beam[i].ray_segment.data.size())
            << "Number of points along beam and entries in P/Q should be the same.";
    }
}

TEST_F(DynamicRaytracingBase, DynamicRayTracingAndKinematicRayTracingShouldResultInSamePath) {
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto ray = rt.trace_ray(initial_state, "TRT");
    auto beam = rt.trace_beam(initial_state, 1, 1, "TRT");
    ASSERT_EQ(ray.size(), beam.size()) << "Size of array and beam are different.";
    for (auto segment_i = 0; segment_i < beam.size(); segment_i++) {
        ASSERT_EQ(beam[segment_i].ray_segment.data.size(), ray[segment_i].data.size())
            << "Different size for ray and beam segment at index " << segment_i;
        for (auto point_i = 0; point_i < ray[segment_i].data.size(); ++point_i) {
            EXPECT_EQ(beam[segment_i].ray_segment.data[point_i], ray[segment_i].data[point_i])
                << "Segments differ at index " << point_i;
        }
    }
}

/**
 * Get absolute path to directory of current source file.
 * @return If the source file is /foo/bar/baz.cpp, return path(/foo/bar)
 */
std::filesystem::path current_source_path() {
    std::filesystem::path p(__FILE__);
    p = p.parent_path();
    return p;
}

template <typename T, typename = void>
struct is_container : std::false_type {};

template <typename T>
struct is_container<
    T, std::conditional_t<
           false,
           std::void_t<typename T::value_type, typename T::size_type, typename T::allocator_type,
                       typename T::iterator, typename T::const_iterator,
                       decltype(std::declval<T>().size()), decltype(std::declval<T>().begin()),
                       decltype(std::declval<T>().end()), decltype(std::declval<T>().cbegin()),
                       decltype(std::declval<T>().cend())>,
           void>> : public std::true_type {};

// override comparison operator in googletest namespace for xtensor shape function that returns
// std::vector or std::array.
namespace testing::internal {
    template <typename Container1, typename Container2,
              std::enable_if_t<::is_container<Container1>::value, int> = 0,
              std::enable_if_t<::is_container<Container2>::value, int> = 0>
    bool operator==(const Container1& c1, const Container2& c2) {
        if (c1.size() != c2.size()) {
            return false;
        }
        for (auto i = 0; i < c1.size(); ++i) {
            if (c1[i] != c2[i]) {
                return false;
            }
        }
        return true;
    }
} // namespace testing::internal


TEST(DynamicRaytracing, TestForRegressionSingleLayer) {
    // trace beam through single layer and compare with previous result.
    VelocityModel model{{{0, 10, 2000, 1}}};
    KinematicRayTracer rt{model};
    auto P_desired =
        xt::load_npy<std::complex<double>>(current_source_path() / "data/P_analytic.npy");
    auto Q_desired =
        xt::load_npy<std::complex<double>>(current_source_path() / "data/Q_analytic.npy");
    // std::cout << P_desired.shape() << std::endl;
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto beam = rt.trace_beam(initial_state, 10, 40, "", 1, 4);
    ASSERT_EQ(beam.segments.size(), 1);
    // result loaded from disk has a different number of steps, compare only first and last entry.
    // compare first/last of P
    EXPECT_EQ(xt::view(beam[0].P, xt::keep(0)), xt::view(xt::squeeze(P_desired), xt::keep(0)));
    EXPECT_EQ(xt::view(beam[0].P, xt::keep(-1)), xt::view(xt::squeeze(P_desired), xt::keep(-1)));
    // compare first/last of Q
    EXPECT_TRUE(xt::allclose(xt::view(beam[0].Q, xt::keep(0)),
                             xt::view(xt::squeeze(Q_desired), xt::keep(0))));
    EXPECT_TRUE(xt::allclose(xt::view(beam[0].Q, xt::keep(-1)),
                             xt::view(xt::squeeze(Q_desired), xt::keep(-1))));
}

#define name_and_value(x) #x << ": " << x << "\n"

TEST(DynamicRayTracing, TestForRegressionMultipleLayers) {
    VelocityModel model{{{0, 10, 2000, 1}, {10, 20, 2000, -1}}};
    KinematicRayTracer rt{model};
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto beam = rt.trace_beam(initial_state, 10, 40, "TRT");
    for (auto i = 0; i < beam.size(); ++i) {
        auto P_desired = xt::load_npy<std::complex<double>>(
            current_source_path() / ("data/P_multilayer" + std::to_string(i) + ".npy"));
        auto Q_desired = xt::load_npy<std::complex<double>>(
            current_source_path() / ("data/Q_multilayer" + std::to_string(i) + ".npy"));
        // compare dimension instead of shape since different number of points may have been
        // calculated during ray tracing.
        EXPECT_EQ(beam[i].P.dimension(), P_desired.dimension());
        EXPECT_EQ(beam[i].Q.dimension(), Q_desired.dimension());
        // compare first and last entry of P and Q (value at interface)
        EXPECT_TRUE(
            xt::allclose(xt::view(beam[i].P, xt::keep(0)), xt::view(P_desired, xt::keep(0))))
            << "First value of P different for beam segment index " << i << ".\n"
            << name_and_value(xt::view(beam[i].P, xt::keep(0)))
            << name_and_value(xt::view(P_desired, xt::keep(0)));
        EXPECT_TRUE(
            xt::allclose(xt::view(beam[i].Q, xt::keep(0)), xt::view(Q_desired, xt::keep(0)), 1E-5, 1E-2))
            << "First value of Q different for beam segment index " << i << ".\n"
            << name_and_value(xt::view(beam[i].Q, xt::keep(0)))
            << name_and_value(xt::view(Q_desired, xt::keep(0)));
        EXPECT_TRUE(
            xt::allclose(xt::view(beam[i].P, xt::keep(-1)), xt::view(P_desired, xt::keep(-1))))
            << "Last value of P different for beam segment index " << i << ".\n"
            << name_and_value(xt::view(beam[i].P, xt::keep(-1)))
            << name_and_value(xt::view(P_desired, xt::keep(-1)));
        EXPECT_TRUE(
            xt::allclose(xt::view(beam[i].Q, xt::keep(-1)), xt::view(Q_desired, xt::keep(-1)), 1E-5, 1E-2))
            << "Last value of Q different for beam segment index " << i << ".\n"
            << name_and_value(xt::view(beam[i].Q, xt::keep(-1)))
            << name_and_value(xt::view(Q_desired, xt::keep(-1)));
    }
}