#ifndef DOUBLEBEAM_CPP_TESTING_UTILS_HPP
#define DOUBLEBEAM_CPP_TESTING_UTILS_HPP

#include "gtest/gtest.h"

/**
 * EXPECT that actual and desired differ by a maximum of atol + rtol * desired.
 * @param actual Actual value
 * @param desired Desired value
 * @param rtol relative tolerance
 * @param atol absolute tolerance
 */
template<typename T>
::testing::AssertionResult Close(T actual, T desired, T rtol = 1E-7, T atol = 0) {
    T max_allowed_difference = atol + rtol * std::abs(desired);
    if (std::abs(actual - desired) <= max_allowed_difference){
        return ::testing::AssertionSuccess();
    } else {
        return ::testing::AssertionFailure() << "Actual value " << actual << " different from desired value " << desired
                                             << " by " << std::abs(actual - desired) << ". Max. allowed difference " << max_allowed_difference;
    }
}

#endif //DOUBLEBEAM_CPP_TESTING_UTILS_HPP
