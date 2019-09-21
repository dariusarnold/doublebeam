#ifndef DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP

#include <array>
#include <functional>

#include "model.hpp"


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
}; // namespace Index


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

/**
 * Wrapper around snells law.
 * @param old_state State at interface before transformation.
 * @param model Velocity model.
 * @param wave_type Which wave type to take ('T' for transmitted, 'R' for reflected).
 * @return New state at interface with changed slowness.
 */
state_type snells_law(const state_type& old_state, const VelocityModel& model, char wave_type);


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
    std::function<double(state_type)>
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
 * Factory function to create a state type for a position and a starting angle.
 * @param x X coordinate of start point of ray.
 * @param y Y  coordinate of start point of ray.
 * @param z Z  coordinate of start point of ray.
 * @param model Model will be evaluated at start point to calculate starting slowness.
 * @param theta Angle against downgoing vertical axis (z) at start point in rad, increasing upwards.
 * Valid range 0 <= theta <= pi.
 * @param phi Angle against x axis at start point in rad, with increasing angle towards the y axis.
 * Valid range 0 <= phi <= 2*pi.
 * @param T Travel time at the start point of the ray.
 * @return state_type with coordinate and slowness calculated from the given parameters.
 */
state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T);


#endif // DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
