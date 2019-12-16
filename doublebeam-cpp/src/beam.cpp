#include "beam.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <utility>


Beam::Beam(Meter beam_width, AngularFrequency beam_frequency, Ray ray, Eigen::Matrix2cd P,
           Eigen::Matrix2cd Q) :
        ray(std::move(ray)),
        segments{BeamSegment(P, Q)},
        m_width(beam_width),
        m_frequency(beam_frequency) {}

Meter Beam::width() const {
    return m_width;
}

AngularFrequency Beam::frequency() const {
    return m_frequency;
}

const BeamSegment& Beam::operator[](int i) const {
    return segments[i];
}

size_t Beam::size() const {
    return segments.size();
}

Slowness Beam::last_slowness() const {
    return ray.last_slowness();
}

Position Beam::last_position() const {
    return ray.last_position();
}

Second Beam::traveltime() const {
    return ray.traveltime();
}

Eigen::Matrix2cd Beam::get_Q(Arclength s) const {
    auto segment_index = find_segment_index(s);
    Velocity v = ray[segment_index].layer_velocity();
    double distance_in_layer = s.length.get() - ray[segment_index].begin().arclength.length.get();
    return segments[segment_index].Q + v.get() * distance_in_layer * segments[segment_index].P;
}

BeamSegment::BeamSegment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q) :
        P(std::move(P)), Q(std::move(Q)) {}

size_t Beam::find_segment_index(Arclength s) const {
    if (s.length.get() < 0) {
        throw std::domain_error(fmt::format("Negative arclength not allowed: {}", s));
    }
    for (size_t i = 0; i < ray.size(); ++i) {
        if (ray[i].begin().arclength.length.get() <= s.length.get() and
            ray[i].end().arclength.length.get() >= s.length.get()) {
            return i;
        }
    }
    throw std::invalid_argument(fmt::format("{} beyond maximum {}.", s, ray.last_arclength()));
}

Arclength Beam::last_arclength() const {
    return ray.last_arclength();
}

const Eigen::Matrix2cd& Beam::last_P() const {
    return segments.back().P;
}

void Beam::add_segment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q) {
    segments.emplace_back(P, Q);
}

Eigen::Matrix2cd Beam::last_Q() const {
    return get_Q(last_arclength());
}
