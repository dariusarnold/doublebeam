#include "ray.hpp"

RaySegment::RaySegment(std::vector<state_type> data, std::vector<double> arclength) :
        data(std::move(data)), arclength(std::move(arclength)) {}


std::vector<RaySegment>::iterator Ray::begin() {
    return segments.begin();
}

std::vector<RaySegment>::iterator Ray::end() {
    return segments.end();
}

std::ptrdiff_t Ray::size() const {
    return segments.size();
}

RaySegment& Ray::operator[](size_t index) {
    return segments[index];
}

const RaySegment& Ray::operator[](size_t index) const {
    return segments[index];
}

slowness_t last_slowness(const Ray& ray) {
    if (ray.size() == 0) {
        throw std::length_error("Accessing empty ray.");
    }
    auto last_segment = ray.segments.back().data.back();
    return {last_segment[Index::PX], last_segment[Index::PY], last_segment[Index::PZ]};
}

position_t last_point(const Ray& ray) {
    if (ray.size() == 0) {
        throw std::length_error("Acessing empty ray.");
    }
    auto last_segment = ray.segments.back().data.back();
    return {last_segment[Index::X], last_segment[Index::Y], last_segment[Index::Z]};
}