#ifndef DOUBLEBEAM_CPP_RAYTRACING_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_HPP

#include "beam.hpp"
#include "model.hpp"
#include "ray.hpp"
#include "raytracing_types.hpp"


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
RayState init_state(Meter x, Meter y, Meter z, const VelocityModel& model, Radian theta, Radian phi,
                    TravelTime T = TravelTime{0_second});

/**
 * Overload for init_state taking tuple of x, y, z coordinate.
 */
RayState init_state(position_t position, const VelocityModel& model, Radian theta, Radian phi,
                    TravelTime T = TravelTime(0_second));

/**
 * Create state type instance.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param z Z coordinate.
 * @param px X component of slowness vector.
 * @param py Y component of slowness vector.
 * @param pz Z component of slowness vector.
 * @param T Travel time.
 * @return
 */
RayState make_state(Meter x, Meter y, Meter z, InverseVelocity px, InverseVelocity py,
                    InverseVelocity pz, TravelTime T = TravelTime(0_second));

/**
 * Overload taking tuples for position and slowness.
 * @param position x, y, z coordinate tuple.
 * @param slowness x, y, z slowness vector.
 * @param T Travel time.
 */
RayState make_state(Position position, Slowness slowness, TravelTime T = TravelTime{0_second});


enum class Status {
    // Ray tracing finished successfully
    Success,
    // Ray left velocity model
    OutOfBounds,
    // Ray reached stop depth
    StopDepthReached
};


template <typename ResultType>
class RayTracingResult {
public:
    Status status;
    // this only holds a ray/beam if ray tracing finished successfully
    std::optional<ResultType> result;

    ResultType value() {
        return result.value();
    }
};


class RayTracer {
public:
    explicit RayTracer(VelocityModel velocity_model);

    /**
     * Trace ray specified by ray code through velocity model.
     * @param initial_state Initial state of the ray.
     * @param ray_code Sequence of characters that decide the ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param stop_depth Depth at which ray tracing is stopped.
     * @return Traced ray.
     */
    RayTracingResult<Ray> trace_ray(const RayState& initial_state,
                                    const std::vector<WaveType>& ray_code = {},
                                    std::optional<Meter> stop_depth = {});

    /**
     * Overload that takes a string and converts it to a ray code.
     * An invalid_argument exception will be thrown when the ray code contains invalid characters.
     */
    RayTracingResult<Ray> trace_ray(const RayState& initial_state, const std::string& ray_code,
                                    std::optional<Meter> stop_depth = {});

    /**
     * Do dynamic ray tracing and return gaussian beam.
     * @param initial_state Initial State (position, slowness, travel time) at the start position.
     * @param beam_width Width of beam.
     * @param beam_frequency Central frequency of beam.
     * @param ray_code Target ray code specifying the ray type to take at an interface.
     * Use 'R' for reflected and 'T' for transmitted.
     * @param stop_depth Depth at which ray tracing is stopped.
     * @return Traced beam.
     */
    RayTracingResult<Beam> trace_beam(const RayState& initial_state, Meter beam_width,
                                      AngularFrequency beam_frequency,
                                      const std::vector<WaveType>& ray_code = {},
                                      std::optional<Meter> stop_depth = {});

    /**
     * Overload that takes a string and converts it to a ray code.
     * An invalid_argument exception will be thrown when the ray code contains invalid characters.
     */
    RayTracingResult<Beam> trace_beam(const RayState& initial_state, Meter beam_width,
                                      AngularFrequency beam_frequency, const std::string& ray_code,
                                      std::optional<Meter> stop_depth = {});

private:
    /**
     * Use analytic ray tracing equation for a constant velocity layer.
     * The equations used are eq. 4 and eq. A.1 from "Analytical ray tracing system: Introducing
     * art, a C-library designed for seismic applications" (Miqueles et al., 2013).
     * @param initial_state State of ray tracing when the ray enters the layer.
     * @param layer Layer which is traced.
     * @return RaySegment for this layer.
     */
    RayTracingResult<RaySegment> trace_layer(const RayState& initial_state, const Layer& layer);

    VelocityModel model;

    std::optional<Meter> stop_depth_m{};
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_HPP
