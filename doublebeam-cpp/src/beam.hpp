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


class BeamSegment {
public:
    BeamSegment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q);

    friend Beam;

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
    Beam(Meter beam_width, AngularFrequency beam_frequency, Ray ray, Eigen::Matrix2cd P,
         Eigen::Matrix2cd Q);

    /**
     * Allow indexing beam to return its segments.
     * @param i
     * @return
     */
    const BeamSegment& operator[](int i) const;

    void add_segment(Eigen::Matrix2cd P, Eigen::Matrix2cd Q);

    [[nodiscard]] size_t size() const;

    [[nodiscard]] Meter width() const;
    [[nodiscard]] AngularFrequency frequency() const;

    [[nodiscard]] Slowness last_slowness() const;

    [[nodiscard]] Position last_position() const;

    [[nodiscard]] Arclength last_arclength() const;

    /**
     * Return value of P at the last beam position.
     */
    [[nodiscard]] const Eigen::Matrix2cd& last_P() const;

    /**
     * Return value of Q at the last beam positon.
     * @return
     */
    [[nodiscard]] Eigen::Matrix2cd last_Q() const;

    [[nodiscard]] Second traveltime() const;

    /**
     * Get matrix Q at specified arclength s.
     * @param s
     * @return
     */
    [[nodiscard]] Eigen::Matrix2cd get_Q(Arclength s) const;

private:
    /**
     * Find index of the segment containing arclength s.
     * @return Index for segment of beam where s(start) <= s <= s(end).
     */
    [[nodiscard]] size_t find_segment_index(Arclength s) const;

    Ray ray;
    std::vector<BeamSegment> segments;
    Meter m_width;
    AngularFrequency m_frequency;
};

#endif // DOUBLEBEAM_CPP_BEAM_HPP
