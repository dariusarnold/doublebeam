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
    auto drt = RayTracer(vm);
    auto initial_state = init_state(0, 0, 0, vm, math::radians(20), 0, 0);
    auto a = std::chrono::high_resolution_clock::now();
    auto beam0 = drt.trace_beam(initial_state, 10, 40, "TTTT");
    auto b = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    std::cout << duration << "\n";
    auto res = measure_runtime([&]() { auto beam = drt.trace_beam(initial_state, 10, 40, "TTTT"); });
    std::cout << res << std::endl;
}