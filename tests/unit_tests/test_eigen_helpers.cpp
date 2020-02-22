#include <gtest/gtest.h>

#include "eigen_helpers.hpp"


class TestTensorAccess : public testing::Test {
protected:
    Eigen::Tensor<int, 3> the_tensor;

    TestTensorAccess() : the_tensor(3, 2, 2) {
        the_tensor.setZero();
        the_tensor(0, 0, 0) = 1;
        the_tensor(0, 1, 0) = 2;
        the_tensor(0, 0, 1) = 3;
        the_tensor(0, 1, 1) = 4;
        the_tensor(2, 0, 0) = -1;
        the_tensor(2, 1, 0) = -2;
        the_tensor(2, 0, 1) = -3;
        the_tensor(2, 1, 1) = -4;
    }
};

TEST_F(TestTensorAccess, TestFirstElement) {
    Eigen::Matrix<int, 2, 2> first_el = first_element(the_tensor);
    Eigen::Matrix<int, 2, 2> expected;
    expected << 1, 2, 3, 4;
    EXPECT_TRUE(first_el == expected) << first_el << "\nnot equal to\n" << expected;
}

TEST_F(TestTensorAccess, TestLastElement) {
    Eigen::Matrix<int, 2, 2> last_el = last_element(the_tensor);
    Eigen::Matrix<int, 2, 2> expected;
    expected << -1, -2, -3, -4;
    EXPECT_TRUE(last_el == expected) << last_el << "\nnot equal to\n" << expected;
}