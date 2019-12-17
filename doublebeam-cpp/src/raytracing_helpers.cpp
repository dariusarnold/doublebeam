#include <algorithm>
#include <cmath>

#include "model.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"


std::tuple<double, double, double> snells_law(double px, double py, double pz, double v_above,
                                              double v_below, WaveType wave_type) {
    double minus_plus = 1.;
    if (wave_type == WaveType::Reflected) {
        return {px, py, -pz};
    } else {
        minus_plus = -1;
    }
    if (not seismo::ray_direction_down(pz)) {
        // for transmitted upgoing ray
        std::swap(v_above, v_below);
    }
    // handle only special case of horizontal interface where normal is vertical.
    // n should be oriented to the side the transmitted wave propagates for the
    // minus_plus relation to work,
    double nz = std::copysign(1., pz);
    // dot product can be simplified since bx and ny are 0, nz is -1
    double p_dot_n = nz * pz;
    double eps = std::copysign(1., p_dot_n);
    // since in the original formula everything subtracted from p is multiplied by n
    // only the pz component changes for horizontal interfaces.
    pz -= (p_dot_n +
           minus_plus * eps *
               std::sqrt(1. / (v_below * v_below) - 1. / (v_above * v_above) + p_dot_n * p_dot_n)) *
          nz;
    return {px, py, pz};
}

Slowness snells_law(const RayState& old_state, const VelocityModel& model, WaveType wave_type) {
    auto [position, slowness, travel_time, arclength] = old_state;
    auto [v_above, v_below] = model.interface_velocities(position.z.get());
    auto [px, py, pz] = snells_law(slowness.px.get(), slowness.py.get(), slowness.pz.get(), v_above,
                                   v_below, wave_type);
    return {InverseVelocity(px), InverseVelocity(py), InverseVelocity(pz)};
}