#include "raytracing.hpp"

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include "model.hpp"
#include "utils.hpp"

namespace odeint = boost::numeric::odeint;


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
    pz -= (p_dot_n +
           minus_plus * eps *
               std::sqrt(1. / (v_below * v_below) - 1. / (v_above * v_above) + p_dot_n * p_dot_n)) *
          nz;
    return {px, py, pz};
}

struct InterfaceCrossed {
    /*
     * Define interface crossing as zero crossing where the function returns
     * values above zero if depth above the interface and below zero if depths
     * below the interface.
     */
    double interface_depth;

    explicit InterfaceCrossed(double interface_depth) : interface_depth(interface_depth){};

    double operator()(const state_type& state) const {
        return interface_depth - state[Index::Z];
    }
};

InterfaceCrossed get_interface_zero_crossing(double pz, const Layer& layer) {
    auto interface_depth = seismo::ray_direction_down(pz) ? layer.bot_depth : layer.top_depth;
    return InterfaceCrossed{interface_depth};
}

/**
 * Integrate a system of ODEs until the Condition has a zero crossing.
 * @tparam System Type of odeint System, callable with signature  sys(x, dxdt, t)-
 * @tparam Condition Type of callable taking a state_type and returning a double.
 * @param x0 Initial state of the system.
 * @param sys System to be integrated.
 * @param cond Return continuous function with a zero crossing for increasing time.
 * @param s_start Start arclength of integration.
 * @param ds Step size of integration.
 * @param max_ds Maximum time step size.
 * @return
 */
template <typename System, typename Condition>
std::pair<std::vector<double>, std::vector<state_type>>
trace_layer(state_type& x0, System sys, Condition cond, double s_start, double ds,
            double max_ds = 1.1) {
    auto stepper =
        odeint::make_dense_output(1.E-10, 1.E-10, max_ds, odeint::runge_kutta_dopri5<state_type>());
    stepper.initialize(x0, s_start, ds);
    std::vector<double> arclengths;
    std::vector<state_type> states;
    // advance stepper until first sign change occurs
    double current_cond = cond(stepper.current_state());
    do {
        states.emplace_back(stepper.current_state());
        arclengths.emplace_back(stepper.current_time());
        stepper.do_step(sys);
    } while (math::same_sign(cond(stepper.current_state()), current_cond));

    // find exact point of zero crossing
    double t0 = stepper.previous_time();
    double t1 = stepper.current_time();
    state_type x_middle;
    boost::uintmax_t max_calls = 1000;
    auto [t_left, t_right] = boost::math::tools::toms748_solve(
        [&](double t) {
            stepper.calc_state(t, x_middle);
            return (cond(x_middle));
        },
        t0, t1, boost::math::tools::eps_tolerance<double>(), max_calls);
    // calculate final position of crossing and save state at crossing
    auto t_middle = (t_left + t_right) * 0.5;
    stepper.calc_state(t_middle, x_middle);
    arclengths.emplace_back(t_middle);
    states.emplace_back(x_middle);
    // workaround for numerical issues where layer boundaries are overshot
    // do this by clamping the value to the interface depth
    auto& last_z = states.back()[Index::Z];
    if (seismo::ray_direction_down(states.back()[Index::PZ])) {
        last_z = std::min(cond.interface_depth, last_z);
    } else {
        last_z = std::max(cond.interface_depth, last_z);
    }
    arclengths.shrink_to_fit();
    states.shrink_to_fit();
    return {arclengths, states};
}


state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T) {
    double velocity = model.eval_at(z);
    auto [px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}

RaySegment::RaySegment(const std::vector<state_type>& states,
                       const std::vector<double>& arclengths) :
        data(states),
        arclength(arclengths) {}


Ray::Ray() : segments(std::vector<RaySegment>()) {}

Ray::Ray(const std::vector<state_type>& states, const std::vector<double>& arclengths) :
        segments{RaySegment{states, arclengths}} {}

Ray::Ray(const RaySegment& segment) : segments(std::vector{segment}) {}


KinematicRayTracer::KinematicRayTracer(VelocityModel velocity_model) :
        model(std::move(velocity_model)) {}

Ray KinematicRayTracer::trace_ray(state_type initial_state, const std::string& ray_code,
                                  double step_size, double max_step) {
    auto layer_index = model.layer_index(initial_state[Index::Z]);
    layer = model[layer_index];
    auto border = get_interface_zero_crossing(initial_state[Index::PZ], layer);
    auto [arclengths, states] = trace_layer(initial_state, *this, border, 0., step_size, max_step);
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
        border = get_interface_zero_crossing(pz_new, layer);
        std::tie(arclengths, states) =
            trace_layer(new_initial_state, *this, border, arclengths.back(), step_size, max_step);
        r.segments.emplace_back(states, arclengths);
    }
    return r;
}

/**
 * Calculate first derivative of inverse of velocity after depth z analytically. 
 * Valid for linear velocity gradient v = v(z) = a * z + b).
 * @param z 
 * @param layer 
 * @return Derivative d/dz of 1/v(z) = -(az+b)^{-2}*a
 */
double dvdz(double z, const Layer& layer) {
    return -layer.gradient /
           ((layer.gradient * z + layer.intercept) * (layer.gradient * z + layer.intercept));
}

void KinematicRayTracer::operator()(const state_type& state, state_type& dxdt,
                                    const double /* s */) {
    auto [x, y, z, px, py, pz, T] = state;
    auto v = model.eval_at(z);
    auto dxds = px * v;
    auto dyds = py * v;
    auto dzds = pz * v;
    auto dpzds = dvdz(z, layer);
    auto dTds = 1. / v;
    dxdt[0] = dxds;
    dxdt[1] = dyds;
    dxdt[2] = dzds;
    dxdt[3] = 0.;
    dxdt[4] = 0.;
    dxdt[5] = dpzds;
    dxdt[6] = dTds;
}
