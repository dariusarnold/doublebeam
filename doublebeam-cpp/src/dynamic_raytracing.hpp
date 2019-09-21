#ifndef DOUBLEBEAM_CPP_DYNAMIC_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_DYNAMIC_RAYTRACING_HPP

#include "beam.hpp"
#include "kinematic_raytracing.hpp"
#include "model.hpp"
#include "raytracing_helpers.hpp"


class DynamicRayTracer {
public:
    explicit DynamicRayTracer(const VelocityModel& model);

    /**
     *
     * @param initial_state
     * @param beam_width
     * @param beam_frequency
     * @param ray_code
     * @param step_size
     * @param max_step
     * @return
     */
    Beam trace_beam(state_type initial_state, double beam_width, double beam_frequency,
                    const std::string& ray_code = "", double step_size = 1., double max_step = 1.1);

private:
    std::pair<xt::xtensor<std::complex<double>, 3>, xt::xtensor<std::complex<double>, 2>>
    transform_PQ_interface(xt::xtensor<std::complex<double>, 2> P,
                           xt::xtensor<std::complex<double>, 2> Q);

    KinematicRayTracer kinematic;
    VelocityModel model;
};


#endif // DOUBLEBEAM_CPP_DYNAMIC_RAYTRACING_HPP
