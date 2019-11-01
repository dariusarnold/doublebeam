#ifndef DOUBLEBEAM_CPP_BEAM_HPP
#define DOUBLEBEAM_CPP_BEAM_HPP

#include <complex>
#include <vector>

#include <xtensor/xtensor.hpp>

#include "ray.hpp"
#include "raytracing_types.hpp"


class BeamSegment {
public:
    BeamSegment(RaySegment segment, xt::xtensor<complex, 3> P, xt::xtensor<complex, 3> Q,
                std::vector<double> v) :
            ray_segment(segment),
            P(P),
            Q(Q),
            v(v) {}
    // Return data of ray segment
    std::vector<state_type> data() const;
    // return arclength of ray segment
    std::vector<double> arclength() const;

    RaySegment ray_segment;
    xt::xtensor<complex, 3> P;
    xt::xtensor<complex, 3> Q;
    std::vector<double> v;
};

/**
 * Get first point (start point) of beam segment.
 * @return x, y, z coordinate triple
 */
position_t first_point(const BeamSegment& bs);

/**
 * Get last point (end point) of beam segment.
 * @return x, y, z coordinate triple
 */
position_t last_point(const BeamSegment& bs);


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
 * Return starting point of beam.
 * @param beam
 * @return
 */
position_t first_point(const Beam& beam);

/**
 * Return last position of a beam.
 * @param beam
 * @return
 */
position_t last_point(const Beam& beam);

// TODO maybe rename to total traveltime
double last_traveltime(const Beam& beam);

double last_arclength(const Beam& beam);

#endif // DOUBLEBEAM_CPP_BEAM_HPP
