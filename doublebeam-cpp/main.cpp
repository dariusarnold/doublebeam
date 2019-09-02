#include <iomanip>
#include "model.h"
#include "printing.h"
#include "integration.hpp"


std::ostream& operator<<(std::ostream& os, const RaySegment& segment){
    os << "Raysegment: " << segment.data.front() << " to " << segment.data.back() << " ( " << segment.data.size() << ") steps";
    return os;
}


int main() {
    std::vector<Layer> layers{{0,   100, 1000, 1},
                              {100, 200, 1100, -1},
                              {200, 300, 1000, 0.5},
                              {300, 400, 1100, 0},
                              {400, 500, 1200, -1}};
    auto vm = VelocityModel(layers);
    auto initial_state = init_state(0., 0., 0., vm, math::radians(20.), math::radians(0.));
    std::vector<state_type> values;
    std::vector<double> times;
    std::cout << std::setprecision(16);
    auto krt = KinematicRayTracer(vm);
    //auto res = measure_runtime([&](){krt.trace_ray(initial_state, "TRT", 1., 1.1);});
    //std::cout << res << "\n";
    auto ray = krt.trace_ray(initial_state, "TTTT", 1., 1.1);
    std::cout << ray.segments << "\n";
    return 0;
}