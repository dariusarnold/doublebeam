#include "doublebeam.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "timing/timing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"


std::ostream& operator<<(std::ostream& os, const RaySegment& segment) {
    os << "Raysegment: " << segment.data.front() << " to " << segment.data.back() << " ( "
       << segment.data.size() << ") steps";
    return os;
}


std::ostream& operator<<(std::ostream& os, const Ray& ray) {
    auto i = 0;
    for (auto& segment : ray.segments) {
        std::cout << i << " " << segment << "\n";
        ++i;
    }
    return os;
}


int main() {
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto db = DoubleBeam(vm);
    auto sources = grid_coordinates(10, 100, 10, 100, 0, 10, 10);
    auto targets = grid_coordinates(10, 100, 10, 100, 450, 10, 10);
    FractureParameters fractures(400, 1, 0, 61, 40, 120, 41);
    db.algorithm(sources, targets, fractures, 10, 40, 0.006);
}