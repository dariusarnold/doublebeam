#include "model.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"
#include <iomanip>
#include <xtensor/xio.hpp>


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
    auto drt = KinematicRayTracer(vm);
    auto initial_state = init_state(0, 0, 0, vm, math::radians(20), 0, 0);
    auto beam = drt.trace_beam(initial_state, 10, 40, "TTTT");
    return 0;
}