#include "ray.hpp"


std::vector<RaySegment>::iterator Ray::begin() {
    return segments.begin();
}

std::vector<RaySegment>::iterator Ray::end() {
    return segments.end();
}
size_t Ray::size() {
    return segments.size();
}
RaySegment& Ray::operator[](size_t index) {
    return segments[index];
}
const RaySegment& Ray::operator[](size_t index) const {
    return segments[index];
}
