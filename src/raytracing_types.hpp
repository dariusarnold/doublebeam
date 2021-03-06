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
#ifndef DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP

#include <array>
#include <complex>
#include <cstddef>
#include <tuple>

#include <boost/container_hash/hash.hpp>
#include <boost/operators.hpp>

#include "units.hpp"


struct Position : boost::equality_comparable<Position> {
    Position() : x(), y(), z() {}
    Position(Meter xx, Meter yy, Meter zz) : x(xx), y(yy), z(zz) {}

    Meter x;
    Meter y;
    Meter z;

    /**
     * Implement hash function because positions are stored in a hash map during unit testing.
     * @param position
     * @return
     */
    friend std::size_t hash_value(const Position& position);
};

/**
 * Calculate distance between point a and b.
 */
inline Meter distance(const Position& a, const Position& b) {
    return Meter(std::sqrt(std::pow(a.x.get() - b.x.get(), 2) + std::pow(a.y.get() - b.y.get(), 2) +
                           std::pow(a.z.get() - b.z.get(), 2)));
}

/**
 * Position is only equality comparable (!=, ==)
 */
bool operator==(const Position& position1, const Position& position2);

std::ostream& operator<<(std::ostream& os, Position position);


struct Slowness : boost::equality_comparable<Slowness> {
    Slowness(InverseVelocity pxx, InverseVelocity pyy, InverseVelocity pzz) :
            px(pxx), py(pyy), pz(pzz) {}

    InverseVelocity px;
    InverseVelocity py;
    InverseVelocity pz;

    /**
     * Change direction of slowness vector to point in the opposite direction.
     */
    inline void flip_direction() {
        px *= -1;
        py *= -1;
        pz *= -1;
    }
};

std::ostream& operator<<(std::ostream& os, Slowness slowness);

/*
 * Slowness is only equality comparable (!=, ==)
 */
bool operator==(const Slowness& slowness1, const Slowness& slowness2);

struct TravelTime : boost::totally_ordered<TravelTime> {
    explicit TravelTime(Second t) : time(t) {}

    Second time;
};

/**
 * Traveltime is totally ordered (<, >, <= ,>=, ==, !=).
 */
bool operator==(const TravelTime& t1, const TravelTime& t2);
bool operator<(const TravelTime& t1, const TravelTime& t2);

std::ostream& operator<<(std::ostream& os, TravelTime travel_time);

struct Arclength : boost::totally_ordered<Arclength> {

    explicit Arclength(Meter arclength) : length(arclength) {}

    Meter length;
};

/**
 * Arclength is totally ordered (<, >, <= ,>=, ==, !=).
 */
bool operator==(const Arclength& l1, const Arclength& l2);
bool operator<(const Arclength& l1, const Arclength& l2);

std::ostream& operator<<(std::ostream& os, Arclength arclength);


struct RayState {
    Position position;
    Slowness slowness;
    TravelTime travel_time;
    Arclength arclength;
};

bool operator==(const RayState& ray_state1, const RayState& ray_state2);
bool operator!=(const RayState& ray_state1, const RayState& ray_state2);

std::ostream& operator<<(std::ostream& os, RayState ray_state);

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
