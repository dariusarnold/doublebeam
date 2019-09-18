#ifndef DOUBLEBEAM_CPP_RAY_HPP
#define DOUBLEBEAM_CPP_RAY_HPP

#include <vector>
#include "raytracing_helpers.hpp"

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
};

#endif // DOUBLEBEAM_CPP_RAY_HPP
