#ifndef DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP

/*
 * Functions and classes used for implementing ray tracing.
 */

#include <array>
#include <functional>

#include <boost/math/tools/toms748_solve.hpp>

#include "model.hpp"
#include "raytracing.hpp"
#include "raytracing_types.hpp"


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
                                              double v_below, WaveType wave_type);

/**
 * Wrapper around snells law.
 * @param old_state State at interface before transformation.
 * @param model Velocity model.
 * @param wave_type Which wave type to take ('T' for transmitted, 'R' for reflected).
 * @return New state at interface with changed slowness.
 */
state_type snells_law(const state_type& old_state, const VelocityModel& model, WaveType wave_type);


class InterfaceCrossed {
public:
    /**
     * Helper class that acts as an event which is triggered once the layer boundary is crossed
     * during ray tracing. If the event is triggered, integration is stopped and the exact location
     * of the layer boundary is found.
     * @param layer Layer in which the ray is traced.
     */
    explicit InterfaceCrossed(const Layer& layer);
    /**
     * Check if depth of current state is outside of the layer.
     * @param state Current state.
     * @return True if interface was crossed with the given state.
     */
    bool operator()(const state_type& state) const;

    /**
     * Return function that has a zero crossing where the layer is crossed. This function is used to
     * locate to exact depth of the interface.
     * @return
     */
    std::function<double(const state_type&)>
    get_zero_crossing_event_function(const state_type& state) const;

    /**
     * Get closest interface depth (either top or bottom) to a given depth.
     * @param state Current state of integration.
     * @return Return top depth of layer if the given state is closer to the top, else return bottom
     * depth.
     */
    double get_closest_layer_depth(const state_type& state) const;

private:
    double top_depth, bottom_depth;
};

/**
 * Helper struct to stop integration once a certain depth has been reached.
 *
 */
class DepthPredicate {
public:
    /**
     * @param z Trace until depth z is reached. Last point of returned ray will then have the exact
     * depth of z.
     */
    explicit DepthPredicate(double z) : z_stop(z) {}
    /**
     * Check if current state has reached the depth.
     * @param state Current raytracing state.
     * @return True if z_stop is between the depths of the previous call and the current call, i.e.
     * if the ray tracing passed the value of z.
     */
    bool operator()(const state_type& state);

private:
    double z_stop;
    std::optional<double> previous_z{};
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
get_state_at_interface(std::function<double(const state_type&)> crossing_function, Stepper stepper) {
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
#endif // DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
