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
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "raytracing_types.hpp"


std::ostream& operator<<(std::ostream& os, Position position) {
    return os << fmt::format("Position({} m, {} m, {} m)", position.x.get(), position.y.get(),
                             position.z.get());
}

std::ostream& operator<<(std::ostream& os, Slowness slowness) {
    return os << fmt::format("Slowness({} s/m, {} s/m, {} s/m)", slowness.px.get(),
                             slowness.py.get(), slowness.pz.get());
}

std::ostream& operator<<(std::ostream& os, TravelTime travel_time) {
    return os << fmt::format("Traveltime({} s)", travel_time.time.get());
}

std::ostream& operator<<(std::ostream& os, Arclength arclength) {
    return os << fmt::format("Arclength({} m)", arclength.length.get());
}

std::ostream& operator<<(std::ostream& os, RayState ray_state) {
    return os << fmt::format("RayState({}, {}, {}, {})", ray_state.position, ray_state.slowness,
                             ray_state.travel_time, ray_state.arclength);
}

bool operator==(const Position& position1, const Position& position2) {
    return position1.x == position2.x and position1.y == position2.y and position1.z == position2.z;
}

bool operator==(const Arclength& length1, const Arclength& length2) {
    return length1.length == length2.length;
}

bool operator<(const Arclength& length1, const Arclength& length2) {
    return length1.length < length2.length;
}

bool operator==(const Slowness& slowness1, const Slowness& slowness2) {
    return slowness1.px == slowness2.px and slowness1.py == slowness2.py and
           slowness1.pz == slowness2.pz;
}

bool operator==(const TravelTime& t1, const TravelTime& t2) {
    return t1.time == t2.time;
}

bool operator<(const TravelTime& t1, const TravelTime& t2) {
    return t1.time < t2.time;
}

bool operator==(const RayState& ray_state1, const RayState& ray_state2) {
    return ray_state1.position == ray_state2.position and
           ray_state1.slowness == ray_state2.slowness and
           ray_state1.travel_time == ray_state2.travel_time and
           ray_state1.arclength == ray_state2.arclength;
}

bool operator!=(const RayState& ray_state1, const RayState& ray_state2) {
    return !(ray_state1 == ray_state2);
}

std::size_t hash_value(const Position& position) {
    std::size_t seed = 0;
    boost::hash_combine(seed, position.x.get());
    boost::hash_combine(seed, position.y.get());
    boost::hash_combine(seed, position.z.get());
    return seed;
}
