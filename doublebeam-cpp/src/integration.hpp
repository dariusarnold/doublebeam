#ifndef DOUBLEBEAM_CPP_INTEGRATION_HPP
#define DOUBLEBEAM_CPP_INTEGRATION_HPP


#include <vector>
#include <array>

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/math/tools/toms748_solve.hpp>

#include "utils.h"

typedef std::array<double, 7> state_type;

namespace odeint = boost::numeric::odeint;

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
template<typename System, typename Condition>
std::pair<std::vector<double>, std::vector<state_type>>
find_crossing(state_type& x0, System sys, Condition cond, double s_start, double ds, double max_ds = 1.1) {
    auto stepper = odeint::make_dense_output(1.E-6, 1.E-6, max_ds, odeint::runge_kutta_dopri5<state_type>());
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
    auto [t_left, t_right] = boost::math::tools::toms748_solve([&](double t){
        stepper.calc_state(t, x_middle); return (cond(x_middle));}, t0, t1,
                boost::math::tools::eps_tolerance<double>(), max_calls);
    // calculate final position of crossing and save state at crossing
    auto t_middle = (t_left + t_right) * 0.5;
    stepper.calc_state(t_middle, x_middle);
    arclengths.emplace_back(t_middle);
    states.emplace_back(x_middle);
    return {arclengths, states};
}


#endif //DOUBLEBEAM_CPP_INTEGRATION_HPP
