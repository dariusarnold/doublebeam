#ifndef DOUBLEBEAM_CPP_TWOPOINT_HPP
#define DOUBLEBEAM_CPP_TWOPOINT_HPP

#include <tuple>

#include <xtensor/xtensor.hpp>

#include "raytracing_types.hpp"
#include "model.hpp"

// TODO evaluate if these type aliases can be moved into TwoPointRayTracing class
using slowness_t = std::tuple<double, double, double>;


class TwoPointRayTracing {
public:
    explicit TwoPointRayTracing(VelocityModel& velocity_model);

    slowness_t trace(position_t source, position_t receiver, double accuracy = 1E-10);

private:
    VelocityModel model;
    size_t num_layers = 0;
    std::vector<double> gradients;
    std::vector<double> intercepts;
    std::vector<double> interface_depths;
    using xtensor_t = xt::xtensor<double, 1, xt::layout_type::row_major>;
    xtensor_t delta_a;
    xtensor_t delta_omega;
    xtensor_t delta_epsilon;
    xtensor_t mu_tilde;
    xtensor_t epsilon_tilde;
    xtensor_t omega_tilde;
    xtensor_t h_tilde;

    auto X_tilde(double q);
    auto X_tilde_prime(double q);
    auto X_tilde_double_prime(double q);
    double f_tilde(double q, double X);
    double f_tilde_prime(double q);
    double f_tilde_double_prime(double q);
    double next_q(double q, double X);
};


#endif // DOUBLEBEAM_CPP_TWOPOINT_HPP
