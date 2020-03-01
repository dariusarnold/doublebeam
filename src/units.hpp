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
#ifndef DOUBLEBEAM_CPP_UNITS_HPP
#define DOUBLEBEAM_CPP_UNITS_HPP

// define macro to get access to constants in Mac clang
#define _USE_MATH_DEFINES
#include <cmath>

#include "strong_types.hpp"


DEFINE_STRONG_TYPE(Meter, double);
DEFINE_TYPE_LITERAL(Meter, _meter);
DEFINE_TYPE_LITERAL_WITH_FACTOR(Meter, _kilometer, 1000);


DEFINE_STRONG_TYPE(Second, double);
DEFINE_TYPE_LITERAL(Second, _second);


DEFINE_STRONG_TYPE(Velocity, double);
DEFINE_TYPE_LITERAL(Velocity, _meter_per_second);
DEFINE_TYPE_LITERAL_WITH_FACTOR(Velocity, _km_per_second, 1000);


DEFINE_STRONG_TYPE(InverseVelocity, double);
DEFINE_TYPE_LITERAL(InverseVelocity, _second_per_meter);


DEFINE_STRONG_TYPE(Radian, double);
DEFINE_TYPE_LITERAL(Radian, _rad);


DEFINE_STRONG_TYPE(Degree, double);
DEFINE_TYPE_LITERAL(Degree, _deg);


/**
 * Convert degree to radians.
 */
inline Radian radians(Degree degree) {
    return Radian(degree.get() * M_PI / 180.);
}

/**
 * Convert radians to degree.
 */
inline Degree degrees(Radian radians) {
    return Degree(radians.get() * 180. / M_PI);
}


DEFINE_STRONG_TYPE(Frequency, double);
DEFINE_TYPE_LITERAL(Frequency, _hertz);

DEFINE_STRONG_TYPE(AngularFrequency, double);
DEFINE_TYPE_LITERAL(AngularFrequency, _rad_per_sec);
DEFINE_TYPE_LITERAL_WITH_FACTOR(AngularFrequency, _angular_from_hertz, 2 * M_PIl);

/**
 * Convert Hertz value to angular frequency.
 */
inline AngularFrequency hertz_to_angular(Frequency freq) {
    return AngularFrequency(2 * M_PI * freq.get());
}

/*
 * Convert angular frequency to Hertz.
 */
inline Frequency angular_to_hertz(AngularFrequency ang) {
    return Frequency(ang.get() / (2 * M_PI));
}

#endif // DOUBLEBEAM_CPP_UNITS_HPP
