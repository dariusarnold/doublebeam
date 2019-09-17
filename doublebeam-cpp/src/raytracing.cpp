#include "raytracing.hpp"

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include "model.hpp"
#include "ray.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"

namespace odeint = boost::numeric::odeint;


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

RaySegment KinematicRayTracer::trace_layer_gradient(const state_type& initial_state,
                                                    const Layer& layer, double s_start, double ds,
                                                    double max_ds) {
    InterfaceCrossed crossing(layer);
    using stepper_t = odeint::runge_kutta_dopri5<state_type>;
    auto stepper = odeint::make_dense_output(1.E-10, 1.E-10, max_ds, stepper_t());
    stepper.initialize(initial_state, s_start, ds);
    std::vector<double> arclengths;
    std::vector<state_type> states;
    // advance stepper until first event occurs
    do {
        states.emplace_back(stepper.current_state());
        arclengths.emplace_back(stepper.current_time());
        stepper.do_step(*this);
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


RaySegment KinematicRayTracer::trace_layer_const(const state_type& initial_state,
                                                 const Layer& layer, double s_start, double ds) {
    auto [x, y, z, px, py, pz, t] = initial_state;
    auto c = layer.intercept;
    auto z_interface = seismo::ray_direction_down(pz) ? layer.bot_depth : layer.top_depth;
    auto s_end = (z_interface - z) / (c * pz);
    size_t num_steps = std::floor(s_end / ds);
    auto s_step = s_end / num_steps;
    // s has to start at 0 because this calculation is done starting from the current point of
    // the ray.
    auto index = 0;
    // one element more since first step (for s = 0) should also be stored.
    std::vector<state_type> states_vector(num_steps + 1);
    std::vector<double> arclengths(num_steps + 1);
    for (auto& el : states_vector) {
        auto arclength_in_layer = index * s_step;
        arclengths[index] = s_start + arclength_in_layer;
        el[Index::X] = x + arclength_in_layer * c * px;
        el[Index::Y] = y + arclength_in_layer * c * py;
        el[Index::Z] = z + arclength_in_layer * c * pz;
        el[Index::PX] = px;
        el[Index::PY] = py;
        el[Index::PZ] = pz;
        el[Index::T] = t + arclength_in_layer / c;
        ++index;
    }
    return {states_vector, arclengths};
}


KinematicRayTracer::KinematicRayTracer(VelocityModel velocity_model) :
        model(std::move(velocity_model)) {}

Ray KinematicRayTracer::trace_ray(state_type initial_state, const std::string& ray_code,
                                  double step_size, double max_step) {
    auto layer_index = model.layer_index(initial_state[Index::Z]);
    current_layer = model[layer_index];
    Ray ray{{trace_layer(initial_state, current_layer, 0., step_size, max_step)}};
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
            current_layer = model[layer_index];
        }
        auto [v_above, v_below] = model.interface_velocities(z);
        auto [px_new, py_new, pz_new] = snells_law(px, py, pz, v_above, v_below, ray_type);
        state_type new_initial_state{x, y, z, px_new, py_new, pz_new, t};
        auto segment = trace_layer(new_initial_state, current_layer,
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
                                    const double /* s */) const {
    auto [x, y, z, px, py, pz, T] = state;
    auto v = model.eval_at(z);
    auto dxds = px * v;
    auto dyds = py * v;
    auto dzds = pz * v;
    auto dpzds = dvdz(z, current_layer);
    auto dTds = 1. / v;
    dfds[0] = dxds;
    dfds[1] = dyds;
    dfds[2] = dzds;
    dfds[3] = 0.;
    dfds[4] = 0.;
    dfds[5] = dpzds;
    dfds[6] = dTds;
}


RaySegment KinematicRayTracer::trace_layer(const state_type& initial_state, const Layer& layer,
                                           double s_start, double ds, double max_ds) {
    if (layer.gradient == 0) {
        return trace_layer_const(initial_state, layer, s_start, ds);
    } else {
        return trace_layer_gradient(initial_state, layer, s_start, ds, max_ds);
    }
}
