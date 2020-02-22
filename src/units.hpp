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
