#ifndef DOUBLEBEAM_CPP_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HPP

#include "beam.hpp"
#include "model.hpp"
#include "ray.hpp"
#include "raytracing_types.hpp"


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
                      double phi, double T = 0);

/**
 * Overload for init_state taking tuple of x, y, z coordinate.
 */
state_type init_state(position_t position, const VelocityModel& model, double theta, double phi, double T = 0);

/**
 * Create state type instance.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param z Z coordinate.
 * @param px X component of slowness vector.
 * @param py Y component of slowness vector.
 * @param pz Z component of slowness vector.
 * @param T Travel time.
 * @return
 */
state_type make_state(double x, double y, double z, double px, double py, double pz, double T = 0);

/**
 * Overload taking tuples for position and slowness.
 * @param position x, y, z coordinate tuple.
 * @param slowness x, y, z slowness vector.
 * @param T Travel time.
 */
state_type make_state(position_t position, slowness_t slowness, double T = 0);


class RayTracer {
public:
    // TODO change to const ref
    explicit RayTracer(VelocityModel velocity_model);

    /**
     * Trace ray specified by ray code through velocity model.
     * @param initial_state Initial state of the ray.
     * @param ray_code Sequence of characters that decide the ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param step_size Initial step size along the ray.
     * @param max_step Maximum step size along the ray.
     * @return Traced ray.
     */
    // TODO take initial state by const ref for ray and beam
    Ray trace_ray(state_type initial_state, const std::vector<WaveType>& ray_code = {},
                  double step_size = 1., double max_step = 1.1);

    /**
     * Overload that takes a string and converts it to a ray code.
     * An invalid_argument exception will be thrown when the ray code contains invalid characters.
     */
    Ray trace_ray(state_type initial_state, const std::string& ray_code, double step_size = 1.,
                  double max_step = 1.1);

    /**
     * Do dynamic ray tracing and return gaussian beam.
     * @param initial_state Initial State (position, slowness, travel time) at the start position.
     * @param beam_width Width of beam.
     * @param beam_frequency Central frequency of beam.
     * @param ray_code Target ray code specifying the ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param step_size Initial step size along the ray.
     * @param max_step Maximum step size along the ray.
     * @return Traced beam.
     */
    Beam trace_beam(state_type initial_state, double beam_width, double beam_frequency,
                    const std::vector<WaveType>& ray_code = {}, double step_size = 1.,
                    double max_step = 1.1);

    /**
     * Overload that takes a string and converts it to a ray code.
     * An invalid_argument exception will be thrown when the ray code contains invalid characters.
     */
    Beam trace_beam(state_type initial_state, double beam_width, double beam_frequency,
                    const std::string& ray_code, double step_size = 1., double max_step = 1.1);

private:
    /**
     * Solve ray tracing equation until layer border is hit.
     * @tparam System Type of odeint System, callable with signature  sys(x, dxdt, t).
     * @param initial_state Initial state of the system.
     * @param sys System to be integrated.
     * @param crossings Callable taking a state_type and returning a bool. Should return true if
     * between the current and the previous call an interface was crossed.
     * @param s_start Start arclength of integration.
     * @param ds Step size of integration.
     * @param max_ds Maximum time step size.
     * @return
     */
    RaySegment trace_layer_gradient(const state_type& initial_state, const Layer& layer,
                                    double s_start, double ds, double max_ds);
    /**
     * Use analytic ray tracing equation for a constant velocity layer.
     * The equations used are eq. 4 and eq. A.1 from "Analytical ray tracing system: Introducing
     * art, a C-library designed for seismic applications" (Miqueles et al., 2013).
     * @param initial_state State of ray tracing when the ray enters the layer.
     * @param layer Layer which is traced.
     * @param s_start Initial value of arc length of ray up to the current point.
     * @param ds Step size of arc length.
     * @return RaySegment for this layer.
     */
    RaySegment trace_layer_const(const state_type& initial_state, const Layer& layer,
                                 double s_start, double ds);
    /**
     * Dispatching function that decides to use analytic or numerical ray tracing depending on
     * whether the layer has constant velocity or not.
     * @param initial_state State of ray tracing when the ray enters the layer.
     * @param layer Layer which is traced.
     * @param s_start Initial value of arc length of ray up to the current point.
     * @param ds Step size of arc length.
     * @param max_ds Maximum time step size.
     * @return RaySegment for this layer.
     */
    RaySegment trace_layer(const state_type& initial_state, const Layer& layer, double s_start,
                           double ds, double max_ds);
    VelocityModel model;
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_HPP
