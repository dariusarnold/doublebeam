/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
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