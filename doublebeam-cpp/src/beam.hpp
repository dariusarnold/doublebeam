#ifndef DOUBLEBEAM_CPP_BEAM_HPP
#define DOUBLEBEAM_CPP_BEAM_HPP

#include <complex>

#include <xtensor/xtensor.hpp>

#include "ray.hpp"


class BeamSegment {
public:
    BeamSegment(RaySegment segment, xt::xtensor<std::complex<double>, 3> P,
                xt::xtensor<std::complex<double>, 3> Q) :
            ray_segment(segment),
            P(P),
            Q(Q) {}
    RaySegment ray_segment;
    xt::xtensor<std::complex<double>, 3> P;
    xt::xtensor<std::complex<double>, 3> Q;
};


class Beam {
public:
    Beam(double beam_width, double beam_frequency);

    std::vector<BeamSegment>::iterator begin();
    std::vector<BeamSegment>::iterator end();
    std::vector<BeamSegment> segments;

private:
    double width;
    double frequency;
};

#endif // DOUBLEBEAM_CPP_BEAM_HPP
