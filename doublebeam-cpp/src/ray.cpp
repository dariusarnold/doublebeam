#include "ray.hpp"


std::size_t Ray::size() const {
    return segments.size();
}

const RaySegment& Ray::operator[](size_t index) const {
    return segments[index];
}

Ray::Ray(RaySegment ray_segment) : segments({ray_segment}) {}

Position Ray::last_position() const {
    if (size() == 0) {
        throw std::length_error("Accessing empty ray.");
    }
    return segments.back().end().position;
}

Slowness Ray::last_slowness() const {
    if (size() == 0) {
        throw std::length_error("Accessing empty ray.");
    }
    return segments.back().end().slowness;
}

void Ray::add_segment(RaySegment ray_segment) {
    segments.push_back(ray_segment);
}

Second Ray::traveltime() const {
    return segments.back().end().travel_time.time;
}

Ray::Ray() : segments() {}

const RayState& Ray::last_state() const {
    return segments.back().end();
}

std::vector<RaySegment>::const_iterator Ray::begin() const {
    return segments.begin();
}

std::vector<RaySegment>::const_iterator Ray::end() const {
    return segments.end();
}

Arclength Ray::last_arclength() const {
    return segments.back().end().arclength;
}

RaySegment::RaySegment(RayState begin, RayState end, Velocity v) :
        begin_m(begin), end_m(end), v_m(v) {}

const RayState& RaySegment::begin() const {
    return begin_m;
}

const RayState& RaySegment::end() const {
    return end_m;
}

Velocity RaySegment::layer_velocity() const {
    return v_m;
}
