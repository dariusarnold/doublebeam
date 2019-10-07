#include "xtensor/xview.hpp"

#include "doublebeam.hpp"
#include "printing.hpp"

/**
 * Calculate horizontal slowness for wave scattered from fractures.
 * @param px X component of slowness vector.
 * @param py Y component of slowness vector.
 * @param phi_hat_x X component of unit vector normal to fracture plane.
 * @param phi_hat_y Y component of unit vector normal to fracture plane.
 * @param fracture_spacing Distance between fracture planes.
 * @param frequency Frequency of wave.
 * @return X, Y component of slowness.
 */
std::tuple<double, double> scattered_slowness(double px, double py, double phi_hat_x,
                                              double phi_hat_y, double fracture_spacing,
                                              double frequency) {
    // pass 0 as pz and phi_hat_z since formula only transforms horizontal slownesses and is only
    // valid for vertical fracture planes.
    auto sign = std::copysign(1., math::dot(px, py, 0, phi_hat_x, phi_hat_y, 0));
    auto px_new = px - sign * phi_hat_x / (fracture_spacing * frequency);
    auto py_new = py - sign * phi_hat_y / (fracture_spacing * frequency);
    return {px_new, py_new};
}


std::vector<position_t> grid_coordinates(double x0, double x1, double y0, double y1, double depth,
                                         int num_x, int num_y) {
    double x_stepsize = std::abs(x1 - x0) / (num_x - 1);
    double y_stepsize = std::abs(y1 - y0) / (num_y - 1);
    std::vector<position_t> points;
    for (auto ix = 0; ix < num_x; ++ix) {
        for (auto iy = 0; iy < num_y; ++iy) {
            points.emplace_back(x0 + ix * x_stepsize, y0 + iy * y_stepsize, depth);
        }
    }
    return points;
}


FractureParameters::FractureParameters(double depth, double phi_hat_x, double phi_hat_y,
                                       int num_fracture_orientations, double spacing_min,
                                       double spacing_max, int num_fracture_spacings) :
        depth(depth),
        phi_hat_x(phi_hat_x),
        phi_hat_y(phi_hat_y),
        orientations(math::generate_vector_arc(num_fracture_orientations, phi_hat_x, phi_hat_y)),
        spacings(math::linspace(spacing_min, spacing_max, num_fracture_spacings)) {}

std::vector<WaveType> direct_ray_code(position_t source, position_t receiver,
                                      const VelocityModel& model) {
    auto n = model.number_of_interfaces_between(std::get<Index::Z>(source),
                                                std::get<Index::Z>(receiver));
    return std::vector<WaveType>(n, WaveType::Transmitted);
}

DoubleBeam::DoubleBeam(const VelocityModel& model) : model(model), twopoint(model), tracer(model) {}

void DoubleBeam::algorithm(std::vector<position_t> source_geometry,
                           std::vector<position_t> target_geometry,
                           FractureParameters fracture_info, double beam_width,
                           double beam_frequency, double window_length) {
    for (const auto& target : target_geometry) {
        for (const auto& source_beam_center : source_geometry) {
            auto slowness = twopoint.trace(target, source_beam_center);
            auto initial_state = make_state(target, slowness, 0);
            auto source_beam =
                tracer.trace_beam(initial_state, beam_width, beam_frequency,
                                  direct_ray_code(source_beam_center, target, model));
            auto last_p = last_slowness(source_beam.value());
            for (auto spacing : fracture_info.spacings) {
                for (auto [phi_hat_x, phi_hat_y] : fracture_info.orientations) {
                    auto a = std::chrono::high_resolution_clock::now();
                    auto [px, py] =
                        scattered_slowness(std::get<0>(last_p), std::get<1>(last_p), phi_hat_x,
                                           phi_hat_y, spacing, beam_frequency);
                    // -pz to reflect beam upwards from target
                    slowness_t new_slowness = {px, py, -std::get<2>(slowness)};
                    initial_state = make_state(target, new_slowness);
                    // reuse ray code since beam should pass through the same layers
                    auto receiver_beam =
                        tracer.trace_beam(initial_state, beam_width, beam_frequency,
                                          direct_ray_code(target, source_beam_center, model));
                    auto b = std::chrono::high_resolution_clock::now();
                    auto duration =
                        std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
                    std::cout << duration << "\n";
                    std::cout << window_length;
                }
            }
        }
    }
}