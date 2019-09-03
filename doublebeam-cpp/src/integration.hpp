#ifndef DOUBLEBEAM_CPP_INTEGRATION_HPP
#define DOUBLEBEAM_CPP_INTEGRATION_HPP

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

state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T = 0);

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

#endif // DOUBLEBEAM_CPP_INTEGRATION_HPP
