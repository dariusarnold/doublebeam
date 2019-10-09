#ifndef DOUBLEBEAM_CPP_BEAM_HPP
#define DOUBLEBEAM_CPP_BEAM_HPP

#include <complex>

#include <xtensor/xtensor.hpp>

#include "ray.hpp"
#include "raytracing_types.hpp"


class BeamSegment {
public:
    BeamSegment(RaySegment segment, xt::xtensor<complex, 3> P, xt::xtensor<complex, 3> Q) :
            ray_segment(segment),
            P(P),
            Q(Q) {}
    RaySegment ray_segment;
    xt::xtensor<complex, 3> P;
    xt::xtensor<complex, 3> Q;
};


class Beam {
public:
    Beam(double beam_width, double beam_frequency, BeamSegment segment);
    Beam(double beam_width, double beam_frequency);

    std::vector<BeamSegment>::iterator begin();
    std::vector<BeamSegment>::iterator end();
    std::vector<BeamSegment> segments;

    /**
     * Allow indexing beam to return its segments.
     * @param i
     * @return
     */
    BeamSegment& operator[](int i);
    const BeamSegment& operator[](int i) const;

    size_t size() const;

    double width() const;
    double frequency() const;

private:
    double m_width;
    double m_frequency;
};

/**
 * Return last slowness value of a beam.
 * @param beam
 * @return
 */
slowness_t last_slowness(const Beam& beam);

/**
 * Return last position of a beam.
 * @param beam 
 * @return 
 */
position_t last_point(const Beam& beam);

#endif // DOUBLEBEAM_CPP_BEAM_HPP
