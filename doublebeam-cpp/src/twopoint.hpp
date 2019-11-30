#ifndef DOUBLEBEAM_CPP_TWOPOINT_HPP
#define DOUBLEBEAM_CPP_TWOPOINT_HPP

#include <tuple>
#include <valarray>

#include "model.hpp"
#include "raytracing_types.hpp"


class TwoPointRayTracing {
public:
    /**
     * Velocity model can only contain layer with constant velocity.
     * @param velocity_model
     */
    explicit TwoPointRayTracing(VelocityModel velocity_model);

    /**
     * Calculate slowness in velocity model to trace from source to receiver.
     *
     * This function is not the pure implementation of the algorithm in the paper
     * "A fast and robust two-point ray tracing method in layered media with constant or linearly
     * varying layer velocity" from Fang and Chen, 2019 but a specialization for a constant layer
     * velocity. Furthermore instead of solving iteratively for q by expanding the offset equation
     * as a Taylor series, a Newton iteration is applied to the Taylor series.
     *
     * @param source x, y, z coordinate tuple
     * @param receiver x, y, z coordinate tuple
     * @param accuracy Iteration is stopped when abs(new_value - current_value) < accuracy.
     * @param max_iterations Iteration is stopped when this number of iterations has been reached.
     * @return slowness px, py, pz
     */
    slowness_t trace(position_t source, position_t receiver, double accuracy = 1E-10,
                     int max_iterations = 20);

    using array_t = std::valarray<double>;

private:
    array_t get_z();

    VelocityModel model;
};


#endif // DOUBLEBEAM_CPP_TWOPOINT_HPP
