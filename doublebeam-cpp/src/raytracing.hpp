#ifndef DOUBLEBEAM_CPP_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HPP

#include <array>
#include <vector>

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


class RaySegment {
public:
    std::vector<state_type> data;
    std::vector<double> arclength;
};

class Ray {
public:
    std::vector<RaySegment> segments;
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
                      double phi, double T = 0);

class KinematicRayTracer {
public:
    explicit KinematicRayTracer(VelocityModel velocity_model);

    /**
     * Trace ray specified by ray code through velocity model.
     * @param initial_state Initial state of the ray.
     * @param ray_code Sequence of characters that decide the ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param step_size Initial step size along the ray.
     * @param max_step Maximum step size along the ray.
     * @return Traced ray.
     */
    Ray trace_ray(state_type initial_state, const std::string& ray_code = "", double step_size = 1.,
                  double max_step = 5.);

    /**
     * This call operator implements the system of ODEs required to compute the ray.
     * @param state
     * @param dxdt
     */
    void operator()(const state_type& state, state_type& dxdt, const double /* s */);

private:
    VelocityModel model;
    Layer layer;
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_HPP
