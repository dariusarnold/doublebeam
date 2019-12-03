#ifndef DOUBLEBEAM_CPP_UNITS_HPP
#define DOUBLEBEAM_CPP_UNITS_HPP

#include "strong_types.hpp"


DEFINE_STRONG_TYPE(Meter, double);
DEFINE_TYPE_LITERAL(Meter, _meter);
DEFINE_TYPE_LITERAL_WITH_FACTOR(Meter, _kilometer, 1000);


DEFINE_STRONG_TYPE(Second, double);


DEFINE_STRONG_TYPE(Velocity, double);
DEFINE_TYPE_LITERAL(Velocity, _meter_per_second);
DEFINE_TYPE_LITERAL_WITH_FACTOR(Velocity, _km_per_second, 1000);


DEFINE_STRONG_TYPE(Radian, double);
DEFINE_TYPE_LITERAL(Radian, _rad);


DEFINE_STRONG_TYPE(Degree, double);
DEFINE_TYPE_LITERAL(Degree , _deg);


DEFINE_STRONG_TYPE(Frequency, double);
DEFINE_TYPE_LITERAL(Frequency , _hertz);


DEFINE_STRONG_TYPE(AngularFrequency, double);
DEFINE_TYPE_LITERAL(AngularFrequency, _rad_per_sec);

#endif // DOUBLEBEAM_CPP_UNITS_HPP
