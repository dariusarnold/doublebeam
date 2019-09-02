#include "integration.hpp"

#include "model.h"
#include "utils.h"


/**
 * Apply snells law to calculate new slowness for horizontal interfaces.
 * @param px X component of slowness vector.
 * @param py Y component of slowness vector.
 * @param pz Z component of slowness vector.
 * @param v_above Velocity on the upper side of the interface
 * @param v_below Velocity on the lower side of the interface.
 * @param wave_type Specify if transmitted ('T') or reflected ('R') wave.
 * @return New slowness values px, py, pz.
 */
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
    pz -= p_dot_n +
          minus_plus * eps * std::sqrt(1. / (v_below * v_below) - 1. / (v_above * v_above) + p_dot_n * p_dot_n) * nz;
    return {px, py, pz};
}


namespace Index {
    /*
     * Indices of variables in state_type.
     * X, Y, Z are cartesian coordinates.
     * PX, PY, PZ are components of slowness vector.
     * T is travel time.
     */
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t Z = 2;
    static constexpr size_t PX = 3;
    static constexpr size_t PY = 4;
    static constexpr size_t PZ = 5;
    static constexpr size_t T = 6;
};


struct InterfaceCrossed {
    /*
     * Define interface crossing as zero crossing where the function returns
     * values above zero if depth above the interface and below zero if depths
     * below the interface.
     */
    double interface_depth;

    explicit InterfaceCrossed(double interface_depth) : interface_depth(interface_depth) {};

    double operator()(const state_type& state) const {
        return interface_depth - state[Index::Z];
    }
};


state_type init_state(double x, double y, double z, const VelocityModel& model,
                      double theta, double phi, double T = 0) {
    const double velocity = model.eval_at(z);
    const auto[px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}

class RaySegment {
public:
    RaySegment(const std::vector<state_type>& states, const std::vector<double>& arclengths) :
            data(states), arclength(arclengths) {}

    std::vector<state_type> data;
    std::vector<double> arclength;
};


class Ray {
public:
    Ray() : segments(std::vector<RaySegment>()) {}

    explicit Ray(const std::vector<state_type>& states, const std::vector<double>& arclengths) :
        segments{RaySegment{states, arclengths}} {}

    explicit Ray(const RaySegment& segment) : segments(std::vector{segment}) {}

    std::vector<RaySegment> segments;
};


class KinematicRayTracer {
public:
    explicit KinematicRayTracer(VelocityModel velocity_model) : model(std::move(velocity_model)) {}

    Ray trace_ray(state_type initial_state, const std::string& ray_code = "", double step_size = 1.,
                  double max_step = 5.) {
        auto layer_index = model.layer_index(initial_state[Index::Z]);
        layer = model[layer_index];
        auto border = get_interface_zero_crossing(initial_state[Index::PZ]);
        auto[arclengths, states] = find_crossing(initial_state, *this, border, 0., step_size, max_step);
        if (ray_code.empty()) {
            return Ray{RaySegment{states, arclengths}};
        }
        Ray r{states, arclengths};
        for (auto ray_type : ray_code) {
            auto [x, y, z, px, py, pz, t] = states.back();
            if (ray_type == 'T') {
                // reflected waves stay in the same layer and dont change the index
                layer_index += seismo::ray_direction_down(pz) ? 1 : -1;
                layer = model[layer_index];
            }
            auto [v_above, v_below] = model.interface_velocities(z);
            auto [px_new, py_new, pz_new] = snells_law(px, py, pz, v_above, v_below, ray_type);
            state_type new_initial_state{x, y, z, px_new, py_new, pz_new, t};
            border = get_interface_zero_crossing(pz_new);
            std::tie(arclengths, states) = find_crossing(new_initial_state, *this, border, arclengths.back(), step_size, max_step);
            r.segments.emplace_back(states, arclengths);
        }
        return r;

    }

    void operator()(const state_type& state, state_type& dxdt, const double /* s */) {
        auto[x, y, z, px, py, pz, T] = state;
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

private:
    InterfaceCrossed get_interface_zero_crossing(double pz) {
        auto interface_depth = seismo::ray_direction_down(pz) ? layer.bot_depth : layer.top_depth;
        return InterfaceCrossed{interface_depth};
    }

    VelocityModel model;

    double dvdz(double z) {
        return -1. / ((layer.gradient * z + layer.intercept) * (layer.gradient * z + layer.intercept));
    }

    Layer layer;
};

