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
DEFINE_TYPE_LITERAL_WITH_FACTOR(AngularFrequency, _angular_from_hertz, 2*M_PIl);

// Those conversion functions here are not beautiful, it would be better if they were members.
inline AngularFrequency hertz_to_angular(Frequency freq) {
    return AngularFrequency(2 * M_PI * freq.get());
}
inline Frequency angular_to_hertz(AngularFrequency ang) {
    return Frequency(ang.get() / (2 * M_PI));
}

#endif // DOUBLEBEAM_CPP_UNITS_HPP
