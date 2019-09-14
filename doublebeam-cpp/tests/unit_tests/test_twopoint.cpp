#include "gtest/gtest.h"
#include <testing_utils.hpp>

#include "model.hpp"
#include "twopoint.hpp"
#include "utils.hpp"

#include <algorithm>
#include <iterator>


class TestTwopointRayTracingBase : public ::testing::Test {
protected:
    TestTwopointRayTracingBase() : twopoint(model) {}
    VelocityModel model = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
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
 */
TEST_P(TestTwoPointRayTracingWithData, CompareWithGivenResult) {
    auto [source, receiver, px_target] = GetParam();
    auto [px, py, pz] = twopoint.trace(source, receiver);
    EXPECT_TRUE(AlmostEqual(px, px_target, 6))
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
    EXPECT_DOUBLE_EQ(pz, 0.00055555555555555556);
}