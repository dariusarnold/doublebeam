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
    auto above = init_state(0_meter, 0_meter, Meter(top), model, 0_rad, 0_rad);
    above.position.z.get() -= 1;
    auto below = init_state(0_meter, 0_meter, Meter(bottom), model, 0_rad, 0_rad);
    below.position.z.get() += 1;
    EXPECT_THROW(rt.trace_beam(above, 0_meter, 0_rad_per_sec), std::domain_error);
    EXPECT_THROW(rt.trace_beam(below, 0_meter, 0_rad_per_sec), std::domain_error);
}

TEST_F(DynamicRaytracingBase, TestBasicCredibilityOfResult) {
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, model, radians(20_deg), 0_rad);
    Meter width(1);
    AngularFrequency frequency(2 * 2 * M_PI);
    Beam beam = rt.trace_beam(initial_state, width, frequency, "TTTT").value();
    EXPECT_EQ(beam.size(), 5) << "Not matching expected number of segments.";
    EXPECT_EQ(beam.width(), width) << "Beam width not set correctly.";
    EXPECT_EQ(beam.frequency(), frequency) << "Beam frequency not set correctly.";
    for (auto i = 0; i < beam.size(); ++i) {
        auto [x, y, z] = beam[i].begin().position;
        EXPECT_EQ(model.eval_at(x, y, z).value(), beam[i].layer_velocity())
            << "Stored and evaluated velocity differ for first point of beam.";
    }
}

TEST_F(DynamicRaytracingBase, DynamicRayTracingAndKinematicRayTracingShouldResultInSamePath) {
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, model, radians(20_deg), 0_rad);
    auto ray = rt.trace_ray(initial_state, "TRT").value();
    auto beam = rt.trace_beam(initial_state, 1_meter, 1_rad_per_sec, "TRT").value();
    ASSERT_EQ(ray.size(), beam.size()) << "Size of array and beam are different.";
    for (auto segment_i = 0; segment_i < beam.size(); segment_i++) {
        const auto& beam_segment = beam[segment_i];
        const auto& ray_segment = ray[segment_i];
        EXPECT_EQ(beam_segment.begin(), ray_segment.begin())
            << "Beam state at begin of beam segment with index " << segment_i
            << " not equal to ray state.";
        EXPECT_EQ(beam_segment.end(), ray_segment.end())
            << "Beam state at end of beam segment with index " << segment_i
            << " not equal to ray state.";
        EXPECT_EQ(beam_segment.layer_velocity(), ray_segment.layer_velocity())
            << "Velocity not equal for ray and beam segment with index " << segment_i;
    }
}


/**
 * Test if index of segment is found correctly from arclength.
 */
struct BeamIndexData {
    Arclength s;
    size_t expected_index;
    friend std::ostream& operator<<(std::ostream& os, const BeamIndexData& x) {
        return os << "BeamIndexData(" << x.s << ", index=" << x.expected_index << ")";
    }
};

class BeamIndexTest : public DynamicRaytracingBase,
                      public testing::WithParamInterface<BeamIndexData> {
protected:
    BeamIndexTest() :
            gauss_beam(rt.trace_beam(initial_state, 10_meter, hertz_to_angular(10_hertz), "TTTT")
                           .value()) {}
    RayState initial_state =
        init_state(0_meter, 0_meter, 0_meter, model, radians(5_deg), radians(0_deg));
    // Arclengths of beam segment
    // begin - end
    // 0 - 100.38198375433474
    //  - 201.06411125891393
    //  - 301.86609181830454
    //  - 402.73176829037379
    //  - 503.69773966239507
    Beam gauss_beam;
};

TEST_P(BeamIndexTest, CompareExpectedAndActualIndex) {
    auto [s, index] = GetParam();
    // Since the index finding function is private we test if the velocity is retrieved correctly.
    EXPECT_EQ(gauss_beam.velocity(s), gauss_beam[index].layer_velocity());
}

INSTANTIATE_TEST_SUITE_P(TestArclengthOfZero, BeamIndexTest,
                         testing::Values(BeamIndexData{Arclength(0_meter), 0}));

INSTANTIATE_TEST_SUITE_P(TestArclengthBeyondArclengthAtBeamEnd, BeamIndexTest,
                         testing::Values(BeamIndexData{Arclength(6000_meter), 4}));

INSTANTIATE_TEST_SUITE_P(TestBeamIndexing, BeamIndexTest,
                         testing::Values(BeamIndexData{Arclength(50_meter), 0},
                                         BeamIndexData{Arclength(150_meter), 1},
                                         BeamIndexData{Arclength(250_meter), 2},
                                         BeamIndexData{Arclength(350_meter), 3},
                                         BeamIndexData{Arclength(450_meter), 4}));


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
    VelocityModel model{{{0_meter, 10_meter, 2000_meter_per_second}}, 1000_meter, 1000_meter};
    RayTracer rt{model};
    auto P_desired =
        load_npy<complex, 3>(current_source_path(__FILE__) / "data/P_analytic_squeeze.npy");
    auto Q_desired =
        load_npy<complex, 3>(current_source_path(__FILE__) / "data/Q_analytic_squeeze.npy");
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, model, radians(20_deg), 0_rad);
    auto beam = rt.trace_beam(initial_state, 10_meter, AngularFrequency(40), "").value();
    ASSERT_EQ(beam.size(), 1);
    // result loaded from disk has a different number of steps, compare only first and last entry.
    // compare first/last of P
    const double low_precision = 1E-2;
    const double high_precision = 1E-10;
    ASSERT_PRED2(matrixCompare(high_precision), beam[0].get_P(), first_element(P_desired));
    ASSERT_PRED2(matrixCompare(high_precision), beam[0].get_P(), last_element(P_desired));
    // compare first/last of Q
    ASSERT_PRED2(matrixCompare(low_precision), beam.get_Q(Arclength{0_meter}),
                 first_element(Q_desired));
    ASSERT_PRED2(matrixCompare(low_precision), beam.last_Q(), last_element(Q_desired));
}

#define name_and_value(x) #x << ": " << x << "\n"

TEST(DynamicRayTracing, TestForRegressionMultipleLayers) {
    VelocityModel model{
        {{0_meter, 10_meter, 2000_meter_per_second}, {10_meter, 20_meter, 2000_meter_per_second}},
        1000_meter,
        1000_meter};
    RayTracer rt{model};
    auto initial_state = init_state(0_meter, 0_meter, 0_meter, model, radians(20_deg), 0_rad);
    auto beam = rt.trace_beam(initial_state, 10_meter, AngularFrequency(40), "TRT").value();
    for (auto i = 0; i < beam.size(); ++i) {
        Eigen::Tensor3cd P_desired = load_npy<std::complex<double>, 3>(
            current_source_path(__FILE__) / ("data/P_multilayer" + std::to_string(i) + ".npy"));
        Eigen::Tensor3cd Q_desired = load_npy<std::complex<double>, 3>(
            current_source_path(__FILE__) / ("data/Q_multilayer" + std::to_string(i) + ".npy"));
        // compare first and last entry of P and Q (value at interface)
        const double low_precision = 1E-2;
        const double high_precision = 1E-10;
        EXPECT_PRED2(matrixCompare(high_precision), beam[i].get_P(), first_element(P_desired))
            << "First value of P different for beam segment index " << i << ".\n"
            << name_and_value(beam[i].get_P()) << name_and_value(first_element(P_desired));
        EXPECT_PRED2(matrixCompare(low_precision), beam[i].Q_begin(), first_element(Q_desired))
            << "First value of Q different for beam segment index " << i << ".\n"
            << name_and_value(beam[i].Q_begin()) << name_and_value(first_element(Q_desired));
        EXPECT_PRED2(matrixCompare(high_precision), beam[i].get_P(), last_element(P_desired))
            << "Last value of P different for beam segment index " << i << ".\n"
            << name_and_value(beam[i].get_P()) << name_and_value(last_element(P_desired));
        EXPECT_PRED2(matrixCompare(low_precision), beam[i].Q_end(), last_element(Q_desired))
            << "Last value of Q different for beam segment index " << i << ".\n"
            << name_and_value(beam[i].Q_end()) << name_and_value(last_element(Q_desired));
    }
}

TEST_F(DynamicRaytracingBase, TestIfFailingCaseWorks) {
    auto initial_state = make_state(
        10_meter, 10_meter, 450_meter, InverseVelocity(-3.8270212473354787e-20),
        InverseVelocity(-0.00062500000000000001), InverseVelocity(-0.00055555555555555556));
    EXPECT_EQ(rt.trace_beam(initial_state, 10_meter, hertz_to_angular(45_hertz), "TTTT").status,
              Status::OutOfBounds);
}