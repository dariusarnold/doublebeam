#include "ray.hpp"


std::vector<RaySegment>::iterator Ray::begin() {
    return segments.begin();
}

std::vector<RaySegment>::iterator Ray::end() {
    return segments.end();
}