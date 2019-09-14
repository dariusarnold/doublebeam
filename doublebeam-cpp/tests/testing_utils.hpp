#ifndef DOUBLEBEAM_CPP_TESTING_UTILS_HPP
#define DOUBLEBEAM_CPP_TESTING_UTILS_HPP

#include "utils.hpp"
#include "gtest/gtest.h"

/**
 * EXPECT that actual and desired differ by a maximum of atol + rtol * desired.
 * @param actual Actual value
 * @param desired Desired value
 * @param rtol relative tolerance
 * @param atol absolute tolerance
 */
template <typename T>
::testing::AssertionResult Close(T actual, T desired, T rtol = 1E-7, T atol = 0) {
    T max_allowed_difference = atol + rtol * std::abs(desired);
    if (std::abs(actual - desired) <= max_allowed_difference) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure()
               << "Actual value " << actual << " different from desired value " << desired << " by "
               << std::abs(actual - desired) << ". Max. allowed difference "
               << max_allowed_difference;
    }
}


/**
 * EXPECT that actual and desired are the same up to a certain number of places (default 7).
 * @tparam T Floating point type.
 * @param actual
 * @param desired
 * @param places Both actual and desired are rounded to this number of places after the decimal
 * point and then compared.
 * @return
 */
template <typename T>
::testing::AssertionResult AlmostEqual(T actual, T desired, int places = 7) {
    if (math::round(actual, places) == math::round(desired, places)) {
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure()
               << actual << " and " << desired << " differ by more than " << places << " digits.";
    }
}

#endif // DOUBLEBEAM_CPP_TESTING_UTILS_HPP