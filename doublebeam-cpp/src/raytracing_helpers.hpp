#ifndef DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP

/*
 * Functions and classes used for implementing ray tracing.
 */

#include <array>
#include <functional>

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
Slowness snells_law(const RayState& old_state, const VelocityModel& model, WaveType wave_type);


#endif // DOUBLEBEAM_CPP_RAYTRACING_HELPERS_HPP
