#include "beam.hpp"


std::vector<BeamSegment>::iterator Beam::begin() {
    return segments.begin();
}

std::vector<BeamSegment>::iterator Beam::end() {
    return segments.end();
}

Beam::Beam(Meter beam_width, AngularFrequency beam_frequency, BeamSegment segment) :
        segments{segment},
        m_width(beam_width),
        m_frequency(beam_frequency) {}

Beam::Beam(Meter beam_width, AngularFrequency beam_frequency) :
        segments(),
        m_width(beam_width),
        m_frequency(beam_frequency) {}

Meter Beam::width() const {
    return m_width;
}

AngularFrequency Beam::frequency() const {
    return m_frequency;
}

BeamSegment& Beam::operator[](int i) {
    return segments[i];
}

const BeamSegment& Beam::operator[](int i) const {
    return segments[i];
}

size_t Beam::size() const {
    return segments.size();
}

slowness_t last_slowness(const Beam& beam) {
    if (beam.size() == 0) {
        throw std::length_error("Accessing empty beam.");
    }
    const auto& last_state = beam.segments.back().ray_segment.data.back();
    return {last_state[Index::PX], last_state[Index::PY], last_state[Index::PZ]};
}

position_t last_point(const Beam& beam) {
    if (beam.size() == 0) {
        throw std::length_error("Acessing empty beam.");
    }
    const auto& last_state = beam.segments.back().ray_segment.data.back();
    return {last_state[Index::X], last_state[Index::Y], last_state[Index::Z]};
}

position_t first_point(const BeamSegment& bs) {
    auto [x, y, z, px, py, pz, t] = bs.data().front();
    return {x, y, z};
}

position_t last_point(const BeamSegment& bs) {
    auto [x, y, z, px, py, pz, t] = bs.data().back();
    return {x, y, z};
}

double last_traveltime(const Beam& beam) {
    return beam.segments.back().data().back()[Index::T];
}

position_t first_point(const Beam& beam) {
    if (beam.size() == 0) {
        throw std::length_error("Acessing empty beam.");
    }
    const auto& first_state = beam.segments.front().ray_segment.data.front();
    return {first_state[Index::X], first_state[Index::Y], first_state[Index::Z]};
}

double last_arclength(const Beam& beam) {
    return beam.segments.back().arclength().back();
}

position_t nth_point(const BeamSegment& bs, size_t n) {
    auto [x, y, z, px, py, pz, t] = bs.ray_segment.data[n];
    return {x, y, z};
}

slowness_t first_slowness(const Beam& beam) {
    if (beam.size() == 0) {
        throw std::length_error("Accessing empty beam.");
    }
    const auto& first_state = beam.segments.front().data().front();
    return {first_state[Index::PX], first_state[Index::PY], first_state[Index::PZ]};
}

std::vector<state_type> BeamSegment::data() const {
    return ray_segment.data;
}

std::vector<double> BeamSegment::arclength() const {
    return ray_segment.arclength;
}
