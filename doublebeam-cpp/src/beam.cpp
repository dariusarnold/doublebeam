#include "beam.hpp"


std::vector<BeamSegment>::iterator Beam::begin() {
    return segments.begin();
}

std::vector<BeamSegment>::iterator Beam::end() {
    return segments.end();
}

Beam::Beam(double beam_width, double beam_frequency, BeamSegment segment) :
        segments{segment},
        m_width(beam_width),
        m_frequency(beam_frequency) {}

Beam::Beam(double beam_width, double beam_frequency) :
        m_width(beam_width),
        m_frequency(beam_frequency) {}

double Beam::width() {
    return m_width;
}

double Beam::frequency() {
    return m_frequency;
}

BeamSegment& Beam::operator[](int i) {
    return segments[i];
}

const BeamSegment& Beam::operator[](int i) const {
    return segments[i];
}

size_t Beam::size() {
    return segments.size();
}