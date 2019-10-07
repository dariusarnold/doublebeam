#ifndef DOUBLEBEAM_CPP_TWOPOINT_HPP
#define DOUBLEBEAM_CPP_TWOPOINT_HPP

#include <tuple>
#include <valarray>

#include "raytracing_types.hpp"
#include "model.hpp"


class TwoPointRayTracing {
public:
    explicit TwoPointRayTracing(const VelocityModel& velocity_model);

    slowness_t trace(position_t source, position_t receiver, double accuracy = 1E-10);

    using array_t = std::valarray<double>;
private:
    VelocityModel model;
    size_t num_layers = 0;
    array_t a;
    array_t b;
    array_t z;
};


#endif // DOUBLEBEAM_CPP_TWOPOINT_HPP
