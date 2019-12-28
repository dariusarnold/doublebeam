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


class Beam;


class BeamSegment : public RaySegment {
public:
    BeamSegment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q, RaySegment segment);

    friend Beam;

    [[nodiscard]] const Eigen::Matrix2cd& get_P() const;

    /**
     * Get Q at begin of ray segment.
     */
    [[nodiscard]] const Eigen::Matrix2cd& Q_begin() const;

    /**
     * Calculate Q at end of ray segment.
     */
    [[nodiscard]] Eigen::Matrix2cd Q_end() const;

private:
    /**
     * Matrix P at begin of beam segment.
     */
    Eigen::Matrix2cd P;
    /**
     * Matrix Q at begin of beam segment.
     */
    Eigen::Matrix2cd Q;
};


class Beam {
public:
    Beam(Meter beam_width, AngularFrequency beam_frequency);

    /**
     * Allow indexing beam to return its segments.
     * @param i
     * @return
     */
    const BeamSegment& operator[](size_t i) const;

    /**
     * Get number of segments.
     * @return
     */
    [[nodiscard]] size_t size() const;

    /**
     * Get last state of the beam, ie the state at the last interface crossing/last ray point.
     * @return
     */
    [[nodiscard]] const RayState& last_state() const;

    /**
     * Get last position of the beam, ie. the position of the end point of the ray.
     */
    [[nodiscard]] const Position& last_position() const;

    /**
     * Get last slowness of a beam, ie. the slowness at the end point of the ray.
     */
    [[nodiscard]] const Slowness& last_slowness() const;

    /**
     * Get last arclength of a beam, ie. the arclength at the end point of a ray.
     */
    [[nodiscard]] Arclength last_arclength() const;

    /**
     * Get total traveltime of the beam from the start point to the end point.
     * @return
     */
    [[nodiscard]] Second traveltime() const;

    /**
     * Enable range based for loop over beam segments.
     * @return
     */
    [[nodiscard]] std::vector<BeamSegment>::const_iterator begin() const;
    [[nodiscard]] std::vector<BeamSegment>::const_iterator end() const;

    void add_segment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q, RaySegment segment);

    [[nodiscard]] Meter width() const;
    [[nodiscard]] AngularFrequency frequency() const;

    /**
     * Return value of P at the last beam position.
     */
    [[nodiscard]] const Eigen::Matrix2cd& last_P() const;

    /**
     * Return value of Q at the last beam positon.
     * @return
     */
    [[nodiscard]] Eigen::Matrix2cd last_Q() const;

    /**
     * Get matrix Q at specified arclength s.
     * @param s
     * @return
     */
    [[nodiscard]] Eigen::Matrix2cd get_Q(Arclength s) const;

    /**
     * Get matrix Q at specified arclength s.
     */
    [[nodiscard]] const Eigen::Matrix2cd& get_P(Arclength s) const;

private:
    /**
     * Find index of the first segment starting from the end of the beam which's begin arclength is
     * less than or equal to arclength s.
     * @return Index for segment of beam where s(start) <= s.
     */
    [[nodiscard]] size_t find_segment_index(Arclength s) const;

    std::vector<BeamSegment> segments;
    Meter m_width;
    AngularFrequency m_frequency;
};

#endif // DOUBLEBEAM_CPP_BEAM_HPP
