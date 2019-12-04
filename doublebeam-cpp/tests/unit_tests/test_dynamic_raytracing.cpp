#include <type_traits>

#include "testing_utils.hpp"
#include <gtest/gtest.h>

#include "eigen_helpers.hpp"
#include "model.hpp"
#include "raytracing.hpp"
#include "utils.hpp"


// version of isApprox that can be used in a ASSERT_PRED/EXPECT_PRED from googletest
auto matrixCompare(double precision) {
    return [=](const auto& matrix1, const auto& matrix2) {
        return matrix1.isApprox(matrix2, precision);
    };
};


template <class Derived>
struct is_eigen : public std::is_base_of<Eigen::DenseBase<Derived>, Derived> {};

// overide stream operator for Eigen matrices
template <class Derived, class = typename std::enable_if<is_eigen<Derived>::value>::type>
::std::ostream& operator<<(::std::ostream& o, const Derived& m) {
    o << "\n" << static_cast<const Eigen::DenseBase<Derived>&>(m);
    return o;
}


class DynamicRaytracingBase : public testing::Test {
protected:
    VelocityModel model = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    RayTracer rt{model};
};

TEST_F(DynamicRaytracingBase, TestThrowWhenStartingOutOfModel) {
    auto [top, bottom] = model.get_top_bottom();
    auto above = init_state(0, 0, top, model, 0, 0, 0);
    above[Index::Z] -= 1;
    auto below = init_state(0, 0, bottom, model, 0, 0, 0);
    below[Index::Z] += 1;
    EXPECT_THROW(rt.trace_beam(above, 0_meter, 0_rad_per_sec), std::domain_error);
    EXPECT_THROW(rt.trace_beam(below, 0_meter, 0_rad_per_sec), std::domain_error);
}

TEST_F(DynamicRaytracingBase, TestBasicCredibilityOfResult) {
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    Meter width(1);
    AngularFrequency frequency(2 * 2 * M_PI);
    auto beam = rt.trace_beam(initial_state, width, frequency, "TTTT").value();
    EXPECT_EQ(beam.segments.size(), 5) << "Not matching expected number of segments.";
    EXPECT_EQ(beam.width(), width) << "Beam width not set correctly.";
    EXPECT_EQ(beam.frequency(), frequency) << "Beam frequency not set correctly.";
    for (auto i = 0; i < beam.segments.size(); ++i) {
        EXPECT_EQ(beam[i].P.dimension(0), beam[i].Q.dimension(0)) << "Length of P and Q different.";
        EXPECT_EQ(beam[i].P.dimension(0), beam[i].ray_segment.data.size())
            << "Number of points along beam and entries in P/Q should be the same.";
        EXPECT_EQ(beam[i].v.size(), beam[i].data().size())
            << "Number of velocities not equal to number of points along beam.";
        auto [x, y, z] = first_point(beam[i]);
        EXPECT_EQ(model.eval_at(x, y, z), beam[i].v.front())
            << "Stored and evaluated velocity differ for first point of beam.";
        std::tie(x, y, z) = last_point(beam[i]);
        EXPECT_EQ(model.eval_at(x, y, z), beam[i].v.back())
            << "Stored and evaluated velocity differ for last point of beam.";
    }
}

TEST_F(DynamicRaytracingBase, DynamicRayTracingAndKinematicRayTracingShouldResultInSamePath) {
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto ray = rt.trace_ray(initial_state, "TRT").value();
    auto beam = rt.trace_beam(initial_state, 1_meter, 1_rad_per_sec, "TRT").value();
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
    VelocityModel model{{{0, 10, 2000, 1}}, 1000, 1000};
    RayTracer rt{model};
    auto P_desired =
        load_npy<complex, 3>(current_source_path(__FILE__) / "data/P_analytic_squeeze.npy");
    auto Q_desired =
        load_npy<complex, 3>(current_source_path(__FILE__) / "data/Q_analytic_squeeze.npy");
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto beam = rt.trace_beam(initial_state, 10_meter, AngularFrequency(40), "", {}, 1, 4).value();
    ASSERT_EQ(beam.segments.size(), 1);
    // result loaded from disk has a different number of steps, compare only first and last entry.
    // compare first/last of P
    const double precision = 1E-6;
    ASSERT_PRED2(matrixCompare(precision), first_element(beam[0].P), first_element(P_desired));
    ASSERT_PRED2(matrixCompare(precision), last_element(beam[0].P), last_element(P_desired));
    // compare first/last of Q
    ASSERT_PRED2(matrixCompare(precision), first_element(beam[0].Q), first_element(Q_desired));
    ASSERT_PRED2(matrixCompare(precision), last_element(beam[0].Q), last_element(Q_desired));
}

#define name_and_value(x) #x << ": " << x << "\n"

TEST(DynamicRayTracing, TestForRegressionMultipleLayers) {
    VelocityModel model{{{0, 10, 2000, 1}, {10, 20, 2000, -1}}, 1000, 1000};
    RayTracer rt{model};
    auto initial_state = init_state(0, 0, 0, model, math::radians(20), 0, 0);
    auto beam = rt.trace_beam(initial_state, 10_meter, AngularFrequency(40), "TRT").value();
    for (auto i = 0; i < beam.size(); ++i) {
        Eigen::Tensor3cd P_desired = load_npy<std::complex<double>, 3>(
            current_source_path(__FILE__) / ("data/P_multilayer" + std::to_string(i) + ".npy"));
        Eigen::Tensor3cd Q_desired = load_npy<std::complex<double>, 3>(
            current_source_path(__FILE__) / ("data/Q_multilayer" + std::to_string(i) + ".npy"));
        // compare dimension instead of shape since different number of points may have been
        // calculated during ray tracing.
        // I know there is a NumDimension, but I get a linker error when using it in the EXPECT_EQ
        // macro.
        EXPECT_EQ(beam[i].P.dimensions().size(), P_desired.dimensions().size());
        EXPECT_EQ(beam[i].Q.dimensions().size(), Q_desired.dimensions().size());
        // compare first and last entry of P and Q (value at interface)
        const double precision = 1E-3;
        EXPECT_PRED2(matrixCompare(precision), first_element(beam[i].P), first_element(P_desired))
            << "First value of P different for beam segment index " << i << ".\n"
            << name_and_value(first_element(beam[i].P)) << name_and_value(first_element(P_desired));
        EXPECT_PRED2(matrixCompare(precision), first_element(beam[i].Q), first_element(Q_desired))
            << "First value of Q different for beam segment index " << i << ".\n"
            << name_and_value(first_element(beam[i].Q)) << name_and_value(first_element(Q_desired));
        EXPECT_PRED2(matrixCompare(precision), last_element(beam[i].P), last_element(P_desired))
            << "Last value of P different for beam segment index " << i << ".\n"
            << name_and_value(last_element(beam[i].P)) << name_and_value(last_element(P_desired));
        EXPECT_PRED2(matrixCompare(precision), last_element(beam[i].Q), last_element(Q_desired))
            << "Last value of Q different for beam segment index " << i << ".\n"
            << name_and_value(last_element(beam[i].Q)) << name_and_value(last_element(Q_desired));
    }
}

TEST_F(DynamicRaytracingBase, TestIfFailingCaseWorks) {
    auto initial_state = make_state(10, 10, 450, -3.8270212473354787e-20, -0.00062500000000000001,
                                    -0.00055555555555555556);
    EXPECT_EQ(rt.trace_beam(initial_state, 10_meter, hertz_to_angular(45), "TTTT").status,
              Status::OutOfBounds);
}