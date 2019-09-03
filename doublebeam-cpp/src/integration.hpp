#ifndef DOUBLEBEAM_CPP_INTEGRATION_HPP
#define DOUBLEBEAM_CPP_INTEGRATION_HPP


#include <vector>
#include <array>

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/math/tools/toms748_solve.hpp>

#include "model.h"
#include "utils.h"


typedef std::array<double, 7> state_type;

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
                                              double v_below, char wave_type);


class RaySegment {
public:
    RaySegment(const std::vector<state_type>& states, const std::vector<double>& arclengths);

    std::vector<state_type> data;
    std::vector<double> arclength;
};


class Ray {
public:
    Ray();

    explicit Ray(const std::vector<state_type>& states, const std::vector<double>& arclengths);

    explicit Ray(const RaySegment& segment);

    std::vector<RaySegment> segments;
};


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
    auto stepper = odeint::make_dense_output(1.E-10, 1.E-10, max_ds, odeint::runge_kutta_dopri5<state_type>());
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
    auto [t_left, t_right] = boost::math::tools::toms748_solve([&](double t) {
                                                                  stepper.calc_state(t, x_middle);
                                                                  return (cond(x_middle));
                                                              }, t0, t1,
                                                              boost::math::tools::eps_tolerance<double>(), max_calls);
    // calculate final position of crossing and save state at crossing
    auto t_middle = (t_left + t_right) * 0.5;
    stepper.calc_state(t_middle, x_middle);
    arclengths.emplace_back(t_middle);
    states.emplace_back(x_middle);
    // workaround for numerical issues where layer boundaries are overshot
    // do this by clamping the value to the interface depth
    auto& last_z = states.back()[Index::Z];
    if (seismo::ray_direction_down(states.back()[Index::PZ])){
        last_z =  std::min(cond.interface_depth, last_z);
    } else {
        last_z = std::max(cond.interface_depth, last_z);
    }
    arclengths.shrink_to_fit();
    states.shrink_to_fit();
    return {arclengths, states};
}


struct InterfaceCrossed {
    /*
     * Define interface crossing as zero crossing where the function returns
     * values above zero if depth above the interface and below zero if depths
     * below the interface.
     */
    double interface_depth;

    explicit InterfaceCrossed(double interface_depth);

    double operator()(const state_type& state) const;
};


state_type init_state(double x, double y, double z, const VelocityModel& model,
                      double theta, double phi, double T = 0);


class KinematicRayTracer {
public:
    explicit KinematicRayTracer(VelocityModel velocity_model);

    Ray trace_ray(state_type initial_state, const std::string& ray_code = "", double step_size = 1.,
                  double max_step = 5.);

    void operator()(const state_type& state, state_type& dxdt, const double /* s */);

private:
    InterfaceCrossed get_interface_zero_crossing(double pz);

    VelocityModel model;

    double dvdz(double z);

    Layer layer;
};


#endif //DOUBLEBEAM_CPP_INTEGRATION_HPP
