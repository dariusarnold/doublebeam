#include "integration.hpp"

#include "model.hpp"
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


InterfaceCrossed::InterfaceCrossed(double interface_depth) : interface_depth(interface_depth){};

double InterfaceCrossed::operator()(const state_type& state) const {
    return interface_depth - state[Index::Z];
}


state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T) {
    double velocity = model.eval_at(z);
    auto [px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}

RaySegment::RaySegment(const std::vector<state_type>& states, const std::vector<double>& arclengths)
: data(states), arclength(arclengths) {}


Ray::Ray() : segments(std::vector<RaySegment>()) {}

Ray::Ray(const std::vector<state_type>& states, const std::vector<double>& arclengths)
: segments{RaySegment{states, arclengths}} {}

Ray::Ray(const RaySegment& segment) : segments(std::vector{segment}) {}


KinematicRayTracer::KinematicRayTracer(VelocityModel velocity_model)
: model(std::move(velocity_model)) {}

Ray KinematicRayTracer::trace_ray(state_type initial_state, const std::string& ray_code,
                                  double step_size, double max_step) {
    auto layer_index = model.layer_index(initial_state[Index::Z]);
    layer = model[layer_index];
    auto border = get_interface_zero_crossing(initial_state[Index::PZ]);
    auto [arclengths, states] =
        find_crossing(initial_state, *this, border, 0., step_size, max_step);
    if (ray_code.empty()) {
        return Ray{RaySegment{states, arclengths}};
    }
    Ray r{states, arclengths};
    for (auto ray_type : ray_code) {
        auto [x, y, z, px, py, pz, t] = states.back();

        if (ray_type == 'T') {
            // reflected waves stay in the same layer and doesn't change the index
            layer_index += seismo::ray_direction_down(pz) ? 1 : -1;
            layer = model[layer_index];
        }
        auto [v_above, v_below] = model.interface_velocities(z);
        auto [px_new, py_new, pz_new] = snells_law(px, py, pz, v_above, v_below, ray_type);
        state_type new_initial_state{x, y, z, px_new, py_new, pz_new, t};
        border = get_interface_zero_crossing(pz_new);
        std::tie(arclengths, states) =
            find_crossing(new_initial_state, *this, border, arclengths.back(), step_size, max_step);
        r.segments.emplace_back(states, arclengths);
    }
    return r;
}

void KinematicRayTracer::operator()(const state_type& state, state_type& dxdt,
                                    const double /* s */) {
    auto [x, y, z, px, py, pz, T] = state;
    auto v = model.eval_at(z);
    auto dxds = px * v;
    auto dyds = py * v;
    auto dzds = pz * v;
    auto dpzds = dvdz(z);
    auto dTds = 1. / v;
    dxdt[0] = dxds;
    dxdt[1] = dyds;
    dxdt[2] = dzds;
    dxdt[3] = 0.;
    dxdt[4] = 0.;
    dxdt[5] = dpzds;
    dxdt[6] = dTds;
}

InterfaceCrossed KinematicRayTracer::get_interface_zero_crossing(double pz) {
    auto interface_depth = seismo::ray_direction_down(pz) ? layer.bot_depth : layer.top_depth;
    return InterfaceCrossed{interface_depth};
}

double KinematicRayTracer::dvdz(double z) {
    return -layer.gradient /
           ((layer.gradient * z + layer.intercept) * (layer.gradient * z + layer.intercept));
}
