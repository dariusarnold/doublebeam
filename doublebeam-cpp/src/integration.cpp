#include <boost/numeric/odeint.hpp>

#include "integration.hpp"
#include "model.h"
#include "utils.h"


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


struct Index {
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


struct InterfaceCrossed {
    /*
     * Define interface crossing as zero crossing where the function returns
     * values above zero if depth above the interface and below zero if depths
     * below the interface.
     */
    double interface_depth;

    explicit InterfaceCrossed(double interface_depth) : interface_depth(interface_depth) {};

    double operator()(const state_type& state) const {
        return interface_depth - state[Index::Z];
    }
};


state_type init_state(double x, double y, double z, const VelocityModel& model,
                      double theta, double phi, double T = 0) {
    const double velocity = model.eval_at(z);
    const auto[px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}