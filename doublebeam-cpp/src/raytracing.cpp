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

class InterfaceCrossed {
public:
    /**
     * Helper class that acts as an event which is triggered once the layer boundary is crossed
     * during ray tracing. If the event is triggered, integration is stopped and the exact location
     * of the layer boundary is found.
     * @param layer Layer in which the ray is traced.
     */
    explicit InterfaceCrossed(const Layer& layer) :
            top_depth(layer.top_depth),
            bottom_depth(layer.bot_depth){};

    /**
     * Check if depth of current state is outside of the layer.
     * @param state Current state.
     * @return True if interface was crossed with the given state.
     */
    bool operator()(const state_type& state) {
        return state[Index::Z] > bottom_depth or state[Index::Z] < top_depth;
    }

    /**
     * Return function that has a zero crossing where the layer is crossed. This function is used to
     * locate to exact depth of the interface.
     * @return
     */
    std::function<double(state_type)> get_zero_crossing_event_function(const state_type& state) {
        if (state[Index::Z] < top_depth) {
            // top interface was crossed
            return [&](const state_type& state) { return top_depth - state[Index::Z]; };
        } else {
            // bottom interface was crossed
            return [&](const state_type& state) { return bottom_depth - state[Index::Z]; };
        }
    }

    /**
     * Get closest interface depth (either top or bottom) to a given depth.
     * @param state Current state of integration.
     * @return Return top depth of layer if the given state is closer to the top, else return bottom
     * depth.
     */
    double get_closest_layer_depth(const state_type& state) {
        return std::abs(state[Index::Z] - top_depth) < std::abs(state[Index::Z] - bottom_depth)
                   ? top_depth
                   : bottom_depth;
    }

private:
    double top_depth, bottom_depth;
};

/**
 * Calculate exact state at interface crossing.
 * @tparam Stepper Boost dense output stepper type.
 * @param crossing_function Continuous function of state that has a zero crossing at the interface
 * depth.
 * @param stepper Stepper used for solving the system of ODEs.
 * @return Pair of: state at interface and arclength at interface.
 */
template <typename Stepper>
std::pair<state_type, double>
get_state_at_interface(std::function<double(state_type)> crossing_function, Stepper stepper) {
    // our integration variable is not time t but arclengths s
    double s0 = stepper.previous_time();
    double s1 = stepper.current_time();
    state_type x_middle;
    boost::uintmax_t max_calls = 1000;
    auto [s_left, s_right] = boost::math::tools::toms748_solve(
        [&](double t) {
            stepper.calc_state(t, x_middle);
            return (crossing_function(x_middle));
        },
        s0, s1, boost::math::tools::eps_tolerance<double>(), max_calls);
    // calculate final position of crossing and save state at crossing
    auto s_middle = (s_left + s_right) * 0.5;
    stepper.calc_state(s_middle, x_middle);
    return {x_middle, s_middle};
}

/**
 * Solve ray tracing equation until layer border is hit.
 * @tparam System Type of odeint System, callable with signature  sys(x, dxdt, t).
 * @param initial_state Initial state of the system.
 * @param sys System to be integrated.
 * @param crossings Callable taking a state_type and returning a bool. Should return true if between
 * the current and the previous call an interface was crossed.
 * @param s_start Start arclength of integration.
 * @param ds Step size of integration.
 * @param max_ds Maximum time step size.
 * @return
 */
template <typename System>
RaySegment trace_layer(state_type& initial_state, System sys, InterfaceCrossed crossing,
                       double s_start, double ds, double max_ds = 1.1) {
    using stepper_t = odeint::runge_kutta_dopri5<state_type>;
    auto stepper = odeint::make_dense_output(1.E-10, 1.E-10, max_ds, stepper_t());
    stepper.initialize(initial_state, s_start, ds);
    std::vector<double> arclengths;
    std::vector<state_type> states;
    // advance stepper until first event occurs
    do {
        states.emplace_back(stepper.current_state());
        arclengths.emplace_back(stepper.current_time());
        stepper.do_step(sys);
    } while (not crossing(stepper.current_state()));
    // find exact point of zero crossing
    auto crossing_function = crossing.get_zero_crossing_event_function(stepper.current_state());
    auto [state_at_crossing, s_at_crossing] = get_state_at_interface(crossing_function, stepper);
    arclengths.emplace_back(s_at_crossing);
    states.emplace_back(state_at_crossing);
    // workaround for numerical issues where layer boundaries are overshot
    // do this by clamping the value to the interface depth
    states.back()[Index::Z] = crossing.get_closest_layer_depth(states.back());
    arclengths.shrink_to_fit();
    states.shrink_to_fit();
    return {states, arclengths};
}

state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T) {
    double velocity = model.eval_at(z);
    auto [px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}

KinematicRayTracer::KinematicRayTracer(VelocityModel velocity_model) :
        model(std::move(velocity_model)) {}

Ray KinematicRayTracer::trace_ray(state_type initial_state, const std::string& ray_code,
                                  double step_size, double max_step) {
    auto layer_index = model.layer_index(initial_state[Index::Z]);
    layer = model[layer_index];
    InterfaceCrossed crossing_events(layer);
    Ray ray{{trace_layer(initial_state, *this, crossing_events, 0., step_size, max_step)}};
    if (ray_code.empty()) {
        return ray;
    }
    for (auto ray_type : ray_code) {
        // TODO this line and the trace_layer call below violates the law of Demeter, try to
        //  refactor it by improving Ray class
        auto [x, y, z, px, py, pz, t] = ray.segments.back().data.back();
        // reflected waves stay in the same layer, so the index doesn't change
        if (ray_type == 'T') {
            layer_index += seismo::ray_direction_down(pz) ? 1 : -1;
            layer = model[layer_index];
        }
        auto [v_above, v_below] = model.interface_velocities(z);
        auto [px_new, py_new, pz_new] = snells_law(px, py, pz, v_above, v_below, ray_type);
        state_type new_initial_state{x, y, z, px_new, py_new, pz_new, t};
        crossing_events = InterfaceCrossed(layer);
        auto segment = trace_layer(new_initial_state, *this, crossing_events,
                                   ray.segments.back().arclength.back(), step_size, max_step);
        ray.segments.push_back(segment);
    }
    return ray;
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

void KinematicRayTracer::operator()(const state_type& state, state_type& dfds,
                                    const double /* s */) {
    auto [x, y, z, px, py, pz, T] = state;
    auto v = model.eval_at(z);
    auto dxds = px * v;
    auto dyds = py * v;
    auto dzds = pz * v;
    auto dpzds = dvdz(z, layer);
    auto dTds = 1. / v;
    dfds[0] = dxds;
    dfds[1] = dyds;
    dfds[2] = dzds;
    dfds[3] = 0.;
    dfds[4] = 0.;
    dfds[5] = dpzds;
    dfds[6] = dTds;
}
