#include <utility>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "beam.hpp"
#include "utils.hpp"


Beam::Beam(Meter beam_width, AngularFrequency beam_frequency) :
        segments(), m_width(beam_width), m_frequency(beam_frequency) {}

Meter Beam::width() const {
    return m_width;
}

AngularFrequency Beam::frequency() const {
    return m_frequency;
}

const BeamSegment& Beam::operator[](size_t i) const {
    return segments[i];
}

Eigen::Matrix2cd Beam::get_Q(Arclength s) const {
    auto segment_index = find_segment_index(s);
    Velocity v = segments[segment_index].layer_velocity();
    double distance_in_layer =
        s.length.get() - segments[segment_index].begin().arclength.length.get();
    return segments[segment_index].Q + v.get() * distance_in_layer * segments[segment_index].P;
}

const Eigen::Matrix2cd& Beam::get_P(Arclength s) const {
    return segments[find_segment_index(s)].P;
}


BeamSegment::BeamSegment(Eigen::Matrix2cd P_, Eigen::Matrix2cd Q_, RaySegment segment) :
        RaySegment(segment.begin(), segment.end(), segment.layer_velocity()),
        P(std::move(P_)),
        Q(std::move(Q_)) {}

size_t Beam::find_segment_index(Arclength s) const {
    // search in reverse direction since beam evaluation is at the end of the ray so we most likely
    // need the last segment.
    auto segment = std::find_if(segments.rbegin(), segments.rend(), [=](const BeamSegment& seg) {
        return seg.begin().arclength <= s;
    });
    if (segment != segments.rend()) {
        // minus one because rend points one past the end and the found segment would be at the
        // actual beginning of the vector which always has distance 1 to the rend but index 0.
        return std::distance(segment, segments.rend()) - 1;;
    }
    throw std::invalid_argument(
        fmt::format("{} is invalid for beam with max arclength of {}.", s, last_arclength()));
}

const Eigen::Matrix2cd& Beam::last_P() const {
    return segments.back().P;
}

void Beam::add_segment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q, RaySegment segment) {
    segments.emplace_back(P, Q, segment);
}

Eigen::Matrix2cd Beam::last_Q() const {
    return get_Q(Beam::last_arclength());
}

std::vector<BeamSegment>::const_iterator Beam::begin() const {
    return segments.begin();
}

std::vector<BeamSegment>::const_iterator Beam::end() const {
    return segments.end();
}

const Eigen::Matrix2cd& BeamSegment::get_P() const {
    return P;
}

const Eigen::Matrix2cd& BeamSegment::Q_begin() const {
    return Q;
}

Eigen::Matrix2cd BeamSegment::Q_end() const {
    return Q + layer_velocity().get() * length().get() * P;
}

size_t Beam::size() const {
    return segments.size();
}

const RayState& Beam::last_state() const {
    return segments.back().end();
}

const Position& Beam::last_position() const {
    return last_state().position;
}

const Slowness& Beam::last_slowness() const {
    return last_state().slowness;
}

Arclength Beam::last_arclength() const {
    return last_state().arclength;
}

Second Beam::traveltime() const {
    return last_state().travel_time.time;
}

Velocity Beam::velocity(Arclength s) const {
    auto index = find_segment_index(s);
    return segments[index].layer_velocity();
}