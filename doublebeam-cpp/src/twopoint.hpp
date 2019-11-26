#ifndef DOUBLEBEAM_CPP_TWOPOINT_HPP
#define DOUBLEBEAM_CPP_TWOPOINT_HPP

#include <tuple>
#include <valarray>

#include "raytracing_types.hpp"
#include "model.hpp"


class TwoPointRayTracing {
public:
    /**
     * Velocity model can only contain layer with constant velocity.
     * @param velocity_model
     */
    explicit TwoPointRayTracing(const VelocityModel& velocity_model);

    /**
     * Calculate slowness in velocity model to trace from source to receiver.
     * @param source x, y, z coordinate tuple
     * @param receiver x, y, z coordinate tuple
     * @param accuracy Iteration is stopped when abs(new_value - current_value) < accuracy.
     * @param max_iterations Iteration is stopped when this number of iterations has been reached.
     * @return slowness px, py, pz
     */
    slowness_t trace(position_t source, position_t receiver, double accuracy = 1E-10, int max_iterations=20);

    using array_t = std::valarray<double>;
private:
    VelocityModel model;
    size_t num_layers = 0;
    array_t b;
    array_t z;
};


#endif // DOUBLEBEAM_CPP_TWOPOINT_HPP
