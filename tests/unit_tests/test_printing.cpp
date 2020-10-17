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

#include "printing.hpp"

#include <iostream>
#include <valarray>


TEST(TestValarray, TestIfPrintedCorrectly) {
    std::valarray<double> a{0.5, 1.5, 2.5, 3};
    std::stringstream ss;
    ss << a;
    ASSERT_EQ(ss.str(), "0.5, 1.5, 2.5, 3");
}

TEST(TestValarray, TestPrintingEmptyArray) {
    std::valarray<int> a;
    std::stringstream ss;
    ss << a;
    ASSERT_EQ(ss.str(), "");
}