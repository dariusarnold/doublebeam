#ifndef DOUBLEBEAM_CPP_RAY_HPP
#define DOUBLEBEAM_CPP_RAY_HPP

#include <vector>
#include "raytracing_types.hpp"


class RaySegment {
public:
    std::vector<state_type> data;
    std::vector<double> arclength;
};

class Ray {
public:
    // TODO make this private
    std::vector<RaySegment> segments;

    std::vector<RaySegment>::iterator begin();

    std::vector<RaySegment>::iterator end();

    /**
     * Return number of ray segments.
     */
    size_t size() const;

    /**
     * Return const reference to a ray segment.
     */
    const RaySegment& operator[](size_t index) const;

    /**
     * Return reference to a ray segment.
     */
    RaySegment& operator[](size_t index);
};

slowness_t last_slowness(const Ray& ray);

position_t last_point(const Ray& ray);

#endif // DOUBLEBEAM_CPP_RAY_HPP
