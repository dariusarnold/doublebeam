#ifndef DOUBLEBEAM_CPP_TWOPOINT_HPP
#define DOUBLEBEAM_CPP_TWOPOINT_HPP

#include <tuple>

using slowness_t = std::tuple<double, double, double>;
using position_t = std::tuple<double, double, double>;

#include "model.hpp"
class TwoPointRayTracing {
public:
    explicit TwoPointRayTracing(VelocityModel& velocity_model);

    slowness_t trace(position_t source, position_t receiver, double accuracy = 1E-10);

private:
    VelocityModel model;
};


#endif // DOUBLEBEAM_CPP_TWOPOINT_HPP
