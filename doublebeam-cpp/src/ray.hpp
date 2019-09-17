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
    std::vector<RaySegment> segments;
};

#endif // DOUBLEBEAM_CPP_RAY_HPP
