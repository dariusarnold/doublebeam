#include "beam.hpp"


Beam::Beam(double beam_width, double beam_frequency) :
        width(beam_width),
        frequency(beam_frequency) {}


std::vector<BeamSegment>::iterator Beam::begin() {
    return segments.begin();
}

std::vector<BeamSegment>::iterator Beam::end() {
    return segments.end();
}
