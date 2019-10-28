#include "xtensor/xview.hpp"
#include <xtensor-blas/xlinalg.hpp>

#include <complex>
#include <numeric>

#include "doublebeam.hpp"
#include "fft.hpp"
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

#define USEDEBUG
#ifdef USEDEBUG
#define msg(x) std::cout << #x << ": " << x << std::endl;
#else
#define msg(x)
#endif

using vector_t = std::tuple<double, double, double>;
struct UnitVectors {
    vector_t e1, e2;
};

UnitVectors get_ray_centred_unit_vectors(const Beam& beam) {
    // This is valid for a planar ray, see 4. from 4.1.3 Cerveny2001. We have a planar ray for a
    // velocity model consisting of horizontal layers, with v = v(z), since nothing changes the
    // slowness along the x or y axis.
    // slowness is in ray plane
    auto [px, py, pz] = last_slowness(beam);
    // other vector in ray plane (connects start and end point)
    auto [x0, y0, z0] = first_point(beam);
    auto [x1, y1, z1] = last_point(beam);
    auto e2 = std::apply(math::normalize, math::cross(px, py, pz, x1 - x0, y1 - y0, z1 - z0));
    auto e1 = std::apply(math::normalize, math::cross(px, py, pz, std::get<0>(e2), std::get<1>(e2),
                                                      std::get<2>(e2)));
    return {e1, e2};
}

std::complex<double> gb_amplitude(const Beam& beam) {
    // this calculates the amplitude at the end of the beam, assuming thats the point you want since
    // the beam reached the surface.
    auto v_s = beam.segments.back().v.back();
    auto v_s0 = beam.segments.front().v.front();
    auto det_Q_s0 = xt::linalg::det(xt::squeeze(xt::view(beam.segments.front().Q, xt::keep(0))));
    auto det_Q_s = xt::linalg::det(xt::squeeze(xt::view(beam.segments.back().Q, xt::keep(-1))));
    return std::sqrt(v_s * det_Q_s0 / (v_s0 * det_Q_s));
}

std::complex<double> gb_exp(const Beam& beam, double q1, double q2) {
    xt::xarray<double> q{q1, q2};
    auto M_s = xt::squeeze(xt::view(beam.segments.back().P, xt::keep(-1))) *
               xt::linalg::inv(xt::squeeze(xt::view(beam.segments.back().Q, xt::keep(-1))));
    return std::exp(std::complex<double>{0, 1} * beam.frequency() *
                    (last_traveltime(beam) + 0.5 * xt::eval(xt::transpose(q) * M_s * q)[0]));
}


std::complex<double> eval_gauss_beam(const Beam& beam, double x, double y, double z) {
    auto [e1, e2] = get_ray_centred_unit_vectors(beam);
    auto [e1x, e1y, e1z] = e1;
    auto [e2x, e2y, e2z] = e2;
    auto e3 = math::cross(e1x, e1y, e1z, e2x, e2y, e2z);
    auto [e3x, e3y, e3z] = e3;
    // auto [xx, yy, zz] = last_point(beam);
    // auto s = last_arclength(beam);
    auto transformation_matrix = math::inv(e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z);
    auto [q1, q2, q3] = math::dot(transformation_matrix, std::make_tuple(x, y, z));
    std::cout << "Ray centred coordinates: " << (impl::Formatter(" ") << q1 << q2 << q3)
              << std::endl;
    std::cout << "A: " << gb_amplitude(beam) << std::endl;
    std::cout << "exp: " << gb_exp(beam, q1, q2);
    return gb_amplitude(beam) * gb_exp(beam, q1, q2);
}


DoubleBeamResult::DoubleBeamResult(size_t num_of_fracture_spacings,
                                   size_t num_of_fracture_orientations) :
        data(num_of_fracture_spacings, num_of_fracture_orientations) {}


DoubleBeamResult DoubleBeam::algorithm(std::vector<position_t> source_geometry,
                                       std::vector<position_t> target_geometry, SeismoData data,
                                       FractureParameters fracture_info, double beam_width,
                                       double beam_frequency,
                                       double __attribute__((unused)) window_length) {
    FFT fft;
    DoubleBeamResult result(fracture_info.spacings.size(), fracture_info.orientations.size());
    for (const auto& target : target_geometry) {
        for (const auto& source_beam_center : source_geometry) {
            auto slowness = twopoint.trace(source_beam_center, target);
            // TODO add overload so declaring initial state is not required for ray tracing
            auto initial_state = make_state(source_beam_center, slowness, 0);
            std::cout << initial_state << std::endl;
            auto source_beam =
                tracer.trace_beam(initial_state, beam_width, beam_frequency,
                                  direct_ray_code(source_beam_center, target, model));
            if (source_beam.status == Status::OutOfBounds) {
                std::cout << "Source beam left model" << std::endl;
                break;
            }
            auto last_p = last_slowness(source_beam.value());
            // calculate scattered slownesses for all fracture spacings/orientations
            for (auto spacing_index = 0U; spacing_index < fracture_info.spacings.size();
                 ++spacing_index) {
                for (auto orientations_index = 0U;
                     orientations_index < fracture_info.orientations.size(); ++orientations_index) {
                    auto [phi_hat_x, phi_hat_y] = fracture_info.orientations[orientations_index];
                    auto [px, py] = scattered_slowness(
                        std::get<0>(last_p), std::get<1>(last_p), phi_hat_x, phi_hat_y,
                        fracture_info.spacings[spacing_index], beam_frequency);
                    // trace receiver beam in scattered direction
                    auto a = std::chrono::high_resolution_clock::now();
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
                    std::cout << duration << " ns"
                              << "\n";
                    if (receiver_beam.status == Status::OutOfBounds) {
                        // beam didn't reach surface, skip
                        std::cout << "Receiver beam left model " << std::endl;
                        break;
                    }
                    auto total_traveltime = last_traveltime(source_beam.value()) +
                                            last_traveltime(receiver_beam.value());
                    // iteration over sources and receivers
                    for (const auto& source_pos : data.sources()) {
                        auto source_beam_val = eval_gauss_beam(source_beam.value(), source_pos.x,
                                                               source_pos.y, source_pos.z);
                        for (const auto& rec_pos : data.receivers()) {
                            auto receiver_beam_val = eval_gauss_beam(source_beam.value(), rec_pos.x,
                                                                     rec_pos.y, rec_pos.z);
                            auto seismogram = cut(data(source_pos, rec_pos), data.timesteps(),
                                                  total_traveltime - window_length / 2,
                                                  total_traveltime + window_length / 2);
                            auto seismogram_freq = fft.execute(seismogram.data);
                            if (seismogram_freq.size() > 1) {
                                std::cerr << "Got " << seismogram_freq.size()
                                          << " frequency bins back.\n";
                            }
                            result.data(spacing_index, orientations_index) +=
                                source_beam_val * receiver_beam_val * seismogram_freq[0];
                        }
                    }
                }
            }
        }
    }
    return result;
}
