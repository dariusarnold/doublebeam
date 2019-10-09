#include <gtest/gtest.h>

#include "printing.hpp"

#include <iostream>
#include <valarray>


TEST(TestValarray, TestIfPrintedCorrectly) {
    std::valarray<double> a{0.5, 1.5, 2.5, 3};
    std::stringstream ss;
    ss << a;
    ASSERT_EQ(ss.str(), "0.5, 1.5, 2.5, 3");
}