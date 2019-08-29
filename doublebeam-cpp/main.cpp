#include "src/integration.cpp"
#include "src/printing.h"


int main() {
    std::vector<Layer> layers{{0,   100, 1000, 1},
                              {100, 200, 1100, -1}};
    auto vm = VelocityModel(layers);
    auto initial_state = init_state(0., 0., 0., vm, math::radians(20.), math::radians(0.));
    std::vector<state_type> values;
    std::vector<double> times;
    auto trace = KinematicRayTracingEquation(vm);

    auto interface_crossed = InterfaceCrossed(100.);

    auto res = find_crossing(initial_state, trace, interface_crossed, 0., 1000., 1.);
    std::cout << res << "\n";

    return 0;
}