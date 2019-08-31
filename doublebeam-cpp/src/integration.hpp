#ifndef DOUBLEBEAM_CPP_INTEGRATION_HPP
#define DOUBLEBEAM_CPP_INTEGRATION_HPP


#include <vector>
#include <array>

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>

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
 * @param precision Precision for zero crossing search. The width of the arclength window
 * in which the zero crossing occurs is decreased until it is smaller than precision.
 * @param max_ds Maximum time step size.
 * @return
 */
template<typename System, typename Condition>
std::pair<std::vector<double>, std::vector<state_type>>
find_crossing(state_type& x0, System sys, Condition cond, const double s_start,
              const double ds, const double precision = 1.E-6, double max_ds = 1.1) {
    auto stepper = odeint::make_dense_output(1.E-6, 1.E-6, max_ds, odeint::runge_kutta_dopri5<state_type>());
    stepper.initialize(x0, s_start, ds);

    // advance stepper until first sign change occurs
    std::vector<double> arclengths;
    std::vector<state_type> states;
    double current_cond = cond(stepper.current_state());
    do {
        states.emplace_back(stepper.current_state());
        arclengths.emplace_back(stepper.current_time());
        stepper.do_step(sys);
    } while (math::same_sign(cond(stepper.current_state()), current_cond));

    // find exact point of zero crossing
    double t0 = stepper.previous_time();
    double t1 = stepper.current_time();
    double t_middle;
    state_type x_middle;
    while (std::abs(t1 - t0) > precision) {
        t_middle = (t0 + t1) * 0.5;
        stepper.calc_state(t_middle, x_middle);
        if (cond(x_middle) == 0) {
            // midpoint found
            break;
        } else if (cond(x_middle) > 0) {
            // we are above the interface crossing, condition change
            // lies below the midpoint
            t0 = t_middle;
        } else {
            // below interface crossing, condition change lies above
            // the midpoint
            t1 = t_middle;
        }
    }
    // interval of size precision; take its midpoint as a final guess
    t_middle = (t0 + t1) * 0.5;
    stepper.calc_state(t_middle, x_middle);
    arclengths.emplace_back(t_middle);
    states.emplace_back(x_middle);
    return {arclengths, states};
}


#endif //DOUBLEBEAM_CPP_INTEGRATION_HPP
