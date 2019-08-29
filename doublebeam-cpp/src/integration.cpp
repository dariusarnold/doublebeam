#include <boost/numeric/odeint.hpp>
#include "model.h"
#include "utils.h"

namespace odeint = boost::numeric::odeint;

typedef std::array<double, 7> state_type;

class KinematicRayTracingEquation {
    const VelocityModel& vm;

    double dvdz(double z) {
        auto layer = vm[vm.layer_index(z)];
        return -1. / ((layer.gradient * z + layer.intercept) * (layer.gradient * z + layer.intercept));
    }

public:
    explicit KinematicRayTracingEquation(const VelocityModel& model) : vm(model) {}

    void operator()(const state_type& state, state_type& dxdt, const double /* s */) {
        auto[x, y, z, px, py, pz, T] = state;
        auto v = vm.eval_at(z);
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
};

enum class Index : size_t {
    /*
     * Indices of variables in state_type.
     * X, Y, Z are cartesian coordinates.
     * PX, PY, PZ are components of slowness vector.
     * T is travel time.
     */
            X = 0,
    Y = 1,
    Z = 2,
    PX = 3,
    PY = 4,
    PZ = 5,
    T = 6,
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
        return interface_depth - state[2];
    }
};


struct store_solver_state {
    std::vector<state_type>& states;
    std::vector<double>& times;

    store_solver_state(std::vector<state_type>& states, std::vector<double>& times) : states(states), times(times) {}

    void operator()(const state_type& x, double t) {
        states.push_back(x);
        times.push_back(t);
    }
};


state_type init_state(double x, double y, double z, const VelocityModel& model,
                      double theta, double phi, double T = 0) {
    const double velocity = model.eval_at(z);
    const auto[px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
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
 * @param precision Precision for zero crossing search. The width of the arclength window
 * in which the zero crossing occurs is decreased until it is smaller than precision.
 * @param max_ds Maximum time step size.
 * @return
 */
template<typename System, typename Condition>
std::vector<std::pair<double, state_type>>
find_crossing(state_type& x0, System sys, Condition cond, const double s_start,
              const double ds, const double precision = 1.E-6, double max_ds = 1.1) {
    auto stepper = odeint::make_dense_output(1.E-6, 1.E-6, max_ds, odeint::runge_kutta_dopri5<state_type>());
    stepper.initialize(x0, s_start, ds);

    // advance stepper until first sign change occurs
    std::vector<std::pair<double, state_type>> states;
    double current_cond = cond(stepper.current_state());
    do {
        states.emplace_back(stepper.current_time(), stepper.current_state());
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
    states.emplace_back(t_middle, x_middle);
    return states;
}
