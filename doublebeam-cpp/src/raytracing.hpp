#ifndef DOUBLEBEAM_CPP_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HPP

#include <array>
#include <vector>

#include "beam.hpp"
#include "model.hpp"
#include "ray.hpp"


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
    // TODO replace ray code string by classes
    Ray trace_ray(state_type initial_state, const std::string& ray_code = "", double step_size = 1.,
                  double max_step = 1.1);

    /**
     * Do dynamic ray tracing and return gaussian beam.
     * @param initial_state Initial State (position, slowness, travel time) at the start position.
     * @param beam_width Width of beam.
     * @param beam_frequency Central frequency of beam.
     * @param ray_code Target ray code specifyinghe ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param step_size Initial step size along the ray.
     * @param max_step Maximum step size along the ray.
     * @return Traced beam.
     */
    Beam trace_beam(state_type initial_state, double beam_width, double beam_frequency,
                    const std::string& ray_code = "", double step_size = 1., double max_step = 1.1);

    /**
     * This call operator implements the system of ODEs required to compute the ray.
     * The method is not called directly from my code, only by the solver.
     * @param state Current state is input from this.
     * @param dfds Next step is stored here.
     * @param s Current arclength along the ray. The ray tracing system of ODEs does not depend
     * on this parameter.
     */
     // TODO make this private if possible to remove it from the api.
    void operator()(const state_type& state, state_type& dfds, const double /* s */) const;

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
    RaySegment trace_layer(const state_type& initial_state, const Layer& layer, double s_start,
                           double ds, double max_ds);
    VelocityModel model;
    Layer current_layer;
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_HPP
