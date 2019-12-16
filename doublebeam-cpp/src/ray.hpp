#ifndef DOUBLEBEAM_CPP_RAY_HPP
#define DOUBLEBEAM_CPP_RAY_HPP

#include <vector>
#include "raytracing_types.hpp"


/**
 * RaySegment for layered constant velocity model where state is only saved at begin/end
 * and the interpolated between them.
 */
class RaySegment {
public:
    RaySegment(RayState begin, RayState end, Velocity v);

    [[nodiscard]] const RayState& begin() const;

    [[nodiscard]] const RayState& end() const;

    /**
     * Get velocity of the layer the ray segment traveled through.
     * @return
     */
    [[nodiscard]] Velocity layer_velocity() const;

private:
    /**
     * State at the beginning of the ray segment.
     */
    RayState begin_m;
    /**
     * State at the end of the ray segment.
     */
    RayState end_m;
    /**
     * Velocity of the layer in which the ray segment traveled.
     */
    Velocity v_m;
};

class Ray {
public:
    /**
     * Create empty ray.
     */
    Ray();

    /**
     * Create ray with one segment.
     * @param ray_segment
     */
    explicit Ray(RaySegment ray_segment);

    /**
     * Append ray segment to end of ray.
     */
    void add_segment(RaySegment ray_segment);

    /**
     * Get last state of the ray, ie the state at the last interface crossing/last ray point.
     * @return
     */
    [[nodiscard]] const RayState& last_state() const;

    /**
     * Get last position of the ray, ie. the position of the end point of the ray.
     */
    [[nodiscard]] Position last_position() const;

    /**
     * Get last slowness of a ray, ie. the slowness at the end point of the ray.
     */
    [[nodiscard]] Slowness last_slowness() const;

    /**
     * Get last arclength of a ray, ie. the arclength at the end point of a ray.
     */
    [[nodiscard]] Arclength last_arclength() const;

    /**
     * Get total traveltime of the ray from the start point to the end point.
     * @return
     */
    [[nodiscard]] Second traveltime() const;

    /**
     * Return number of ray segments.
     */
    [[nodiscard]] std::size_t size() const;

    /**
     * Enable ranged-for loops over ray segments.
     */
    [[nodiscard]] std::vector<RaySegment>::const_iterator begin() const;
    [[nodiscard]] std::vector<RaySegment>::const_iterator end() const;

    /**
     * Return certain ray segment.
     * @param index Index starting at 0 for the first segment (containing the ray start point),
     * ending at N-1 (containing the ray end point) for a ray with N segments.
     * @return Ray segment with index i.
     */
    const RaySegment& operator[](size_t index) const;

private:
    std::vector<RaySegment> segments;

};

#endif // DOUBLEBEAM_CPP_RAY_HPP
