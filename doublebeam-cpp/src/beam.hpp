#ifndef DOUBLEBEAM_CPP_BEAM_HPP
#define DOUBLEBEAM_CPP_BEAM_HPP

#include <complex>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "eigen_helpers.hpp"
#include "ray.hpp"
#include "raytracing_types.hpp"
#include "units.hpp"


class BeamSegment {
public:
    BeamSegment(RaySegment segment, const Eigen::Tensor3cd& P, const Eigen::Tensor3cd& Q,
                double v) :
            ray_segment(segment), P(std::move(P)), Q(std::move(Q)), v(v) {}
    // Return data of ray segment
    std::vector<state_type> data() const;
    // return arclength of ray segment
    std::vector<double> arclength() const;

    RaySegment ray_segment;
    Eigen::Tensor3cd P;
    Eigen::Tensor3cd Q;
    double v;
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

/**
 * Get nth point of beam segment.
 * @param The index n increases in the beam direction and starts at 0. There is no bounds checking.
 * @return x, y, z coordinate triple
 */
position_t nth_point(const BeamSegment& bs, size_t n);


class Beam {
public:
    Beam(Meter beam_width, AngularFrequency beam_frequency, BeamSegment segment);
    Beam(Meter beam_width, AngularFrequency beam_frequency);

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

    Meter width() const;
    AngularFrequency frequency() const;

private:
    Meter m_width;
    AngularFrequency m_frequency;
};

/**
 * Return last slowness value of a beam.
 * @param beam
 * @return
 */
slowness_t last_slowness(const Beam& beam);

/**
 * Return first slowness value of a beam.
 * @param beam
 * @return
 */
slowness_t first_slowness(const Beam& beam);

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
