#ifndef DOUBLEBEAM_CPP_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HPP

#include <array>
#include <vector>

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include "model.hpp"
#include "utils.hpp"

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

namespace odeint = boost::numeric::odeint;

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

    void operator()(const state_type& state, state_type& dxdt, const double /* s */);

private:
    VelocityModel model;
    Layer layer;
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_HPP
