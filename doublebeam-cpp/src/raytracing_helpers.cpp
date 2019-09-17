#include <algorithm>
#include <cmath>

#include "model.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"


std::tuple<double, double, double> snells_law(double px, double py, double pz, double v_above,
                                              double v_below, char wave_type) {
    double minus_plus = 1.;
    if (wave_type == 'R') {
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


InterfaceCrossed::InterfaceCrossed(const Layer& layer) :
        top_depth(layer.top_depth),
        bottom_depth(layer.bot_depth){};

bool InterfaceCrossed::operator()(const state_type& state) const {
    return state[Index::Z] > bottom_depth or state[Index::Z] < top_depth;
}

std::function<double(state_type)>
InterfaceCrossed::get_zero_crossing_event_function(const state_type& state) const {
    if (state[Index::Z] < top_depth) {
        // top interface was crossed
        return [&](const state_type& state) { return top_depth - state[Index::Z]; };
    } else {
        // bottom interface was crossed
        return [&](const state_type& state) { return bottom_depth - state[Index::Z]; };
    }
}

double InterfaceCrossed::get_closest_layer_depth(const state_type& state) const {
    return std::abs(state[Index::Z] - top_depth) < std::abs(state[Index::Z] - bottom_depth)
               ? top_depth
               : bottom_depth;
}


state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T) {
    double velocity = model.eval_at(z);
    auto [px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}
