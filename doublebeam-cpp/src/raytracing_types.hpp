#ifndef DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP

#include <array>
#include <complex>
#include <cstddef>
#include <tuple>

#include "units.hpp"


struct Position {
    Meter x;
    Meter y;
    Meter z;
};

std::ostream& operator<<(std::ostream& os, Position position);

struct Slowness {
    InverseVelocity px;
    InverseVelocity py;
    InverseVelocity pz;
};

std::ostream& operator<<(std::ostream& os, Slowness slowness);

struct TravelTime {
    Second time;
};

std::ostream& operator<<(std::ostream& os, TravelTime travel_time);

struct Arclength {
    Meter length;
};

std::ostream& operator<<(std::ostream& os, Arclength arclength);

struct RayState {
    Position position;
    Slowness slowness;
    TravelTime travel_time;
    Arclength arclength;
};

std::ostream& operator<<(std::ostream& os, RayState ray_state);

using position_t = std::tuple<double, double, double>;

using slowness_t = std::tuple<double, double, double>;

using complex = std::complex<double>;

/**
 * Specifies wave type to take when crossing an interface. Can be used to build ray codes (sequences
 * of wave types), which specify the behaviour of a ray during it's travel through the velocity
 * model.
 */
enum class WaveType : char { Transmitted = 'T', Reflected = 'R' };

inline WaveType to_wavetype(char c) {
    if (c == static_cast<char>(WaveType::Transmitted)) {
        return WaveType::Transmitted;
    } else if (c == static_cast<char>(WaveType::Reflected)) {
        return WaveType::Reflected;
    } else {
        std::stringstream ss;
        ss << "Wave type specifier " << c
           << " is not valid. Valid are: T (transmitted), R (reflected).";
        throw std::runtime_error(ss.str());
    }
}

inline char to_char(WaveType wave_type) {
    return static_cast<char>(wave_type);
}

#endif // DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP
