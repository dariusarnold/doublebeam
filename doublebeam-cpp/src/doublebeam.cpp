#include <complex>

#include <fmt/format.h>

#include "doublebeam.hpp"
#include "eigen_helpers.hpp"
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
                                              Frequency frequency) {
    // pass 0 as pz and phi_hat_z since formula only transforms horizontal slownesses and is only
    // valid for vertical fracture planes.
    auto sign = std::copysign(1., math::dot(px, py, 0, phi_hat_x, phi_hat_y, 0));
    auto px_new = px - sign * phi_hat_x / (fracture_spacing * frequency.get());
    auto py_new = py - sign * phi_hat_y / (fracture_spacing * frequency.get());
    return {px_new, py_new};
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
    auto n = model.number_of_interfaces_between(std::get<2>(source), std::get<2>(receiver));
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

/**
 * Return orthogonal ray centred unit vectors e1, e2 at the last point of the beam.
 * @param beam
 * @return
 */
UnitVectors get_ray_centred_unit_vectors(const Beam& beam) {
    // This is valid for a planar ray, see 4. from 4.1.3 Cerveny2001. We have a planar ray for a
    // velocity model consisting of horizontal layers, with v = v(z), since nothing changes the
    // slowness along the x or y axis.
    // In this simple case we have to choose one vector (e2) perpendicular to the plane at the
    // initial point of the beam. The vector e2 will then be constant along the ray.
    // The second unit vector perpendicular to the ray can be found by a cross product with e3,
    // which has the same direction as the slowness.
    // slowness is in ray plane
    auto [px, py, pz] = beam.last_slowness();
    if (px.get() == 0 and py.get() == 0) {
        // special case of purely up/down going ray, infinite number of planes contain this ray.
        // Just select a set of unit vectors.
        return {{1, 0, 0}, {0, 1, 0}};
    }
    // Select unit vector perpendicular to ray plane. While Cerveny2001 mentions selecting e2 to be
    // perpendicular at the initial point, we can also select it at the last point (as done here)
    // since e2 will be constant (as said by Cerveny2001 as well).
    auto e2 = math::cross(0, 0, 1, px.get(), py.get(), 0);
    // find other unit vector
    auto e1 = math::cross(-px.get(), -py.get(), -pz.get(), std::get<0>(e2), std::get<1>(e2),
                          std::get<2>(e2));
    e2 = std::apply(math::normalize, e2);
    e1 = std::apply(math::normalize, e1);
    return {e1, e2};
}

std::complex<double> gb_amplitude(const Beam& beam) {
    // this calculates the amplitude at the end of the beam, assuming thats the point you want since
    // the beam reached the surface.
    // Since velocity is constant in one layer its factored out
    std::complex<double> det_Q_s0 = beam.get_Q(Arclength{0_meter}).determinant();
    std::complex<double> det_Q_s = beam.get_Q(beam.last_arclength()).determinant();
    return std::sqrt(det_Q_s0 / det_Q_s);
}

std::complex<double> gb_exp(const Beam& beam, double q1, double q2,
                            std::complex<double>& complex_traveltime) {
    using namespace std::complex_literals;
    Eigen::Vector2d q{q1, q2};
    auto M_s = beam.last_P() * beam.last_Q().inverse();
    complex_traveltime = beam.traveltime().get() + 0.5 * (q.transpose() * M_s * q)[0];
    return std::exp(1i * beam.frequency().get() * complex_traveltime);
}

double unit_vect{0};
double amplt{0};
double expt{0};
double restt{0};

struct BeamEvalResult {
    std::complex<double> gb_value;
    std::complex<double> complex_traveltime;
};

BeamEvalResult eval_gauss_beam(const Beam& beam, double x, double y, double z) {
    using namespace std::complex_literals;
    auto [e1, e2] = get_ray_centred_unit_vectors(beam);
    auto [e1x, e1y, e1z] = e1;
    auto [e2x, e2y, e2z] = e2;
    auto e3 = math::cross(e1x, e1y, e1z, e2x, e2y, e2z);
    e3 = math::scale_vector(e3, 1);
    auto [e3x, e3y, e3z] = e3;
    // We want to convert from general Cartesian system to ray centred, use right eq. 4.1.31 from
    // Cerveny2001.
    auto transformation_matrix = std::make_tuple(e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z);
    auto [q1, q2, q3] = math::dot(transformation_matrix, std::make_tuple(x, y, z));
    auto amp = gb_amplitude(beam);
    std::complex<double> complex_traveltime;
    auto exp = gb_exp(beam, q1, q2, complex_traveltime);
    return {std::conj(amp * exp), complex_traveltime};
}

template <typename T>
BeamEvalResult eval_gauss_beam(const Beam& beam, T cartesian_position) {
    return eval_gauss_beam(beam, cartesian_position.x, cartesian_position.y, cartesian_position.z);
}


DoubleBeamResult::DoubleBeamResult(size_t num_of_fracture_spacings,
                                   size_t num_of_fracture_orientations) :
        data(num_of_fracture_spacings, num_of_fracture_orientations) {
    data.setZero();
}

double cutt = 0.;
double fftt = 0;
double evalt = 0;
double beamt = 0;

double squared_distance(const Position& x, const PositionWithIndex& pos) {
    auto [x0, x1, x2] = x;
    return std::pow(x0.get() - pos.x, 2) + std::pow(x1.get() - pos.y, 2) +
           std::pow(x2.get() - pos.z, 2);
}

std::complex<double> stack(const Beam& source_beam, const Beam& receiver_beam,
                           const SeismoData& data, double window_length, double max_eval_distance) {
    std::complex<double> stacking_result(0, 0);
    std::vector<BeamEvalResult> source_beam_values, receiver_beam_values;
    source_beam_values.reserve(std::size(data.sources()));
    receiver_beam_values.reserve(std::size(data.receivers()));
    using namespace std::placeholders;
    auto a = std::chrono::high_resolution_clock::now();
    std::transform(data.sources().begin(), data.sources().end(),
                   std::back_inserter(source_beam_values),
                   std::bind(eval_gauss_beam<Source>, source_beam, _1));
    std::transform(data.receivers().begin(), data.receivers().end(),
                   std::back_inserter(receiver_beam_values),
                   std::bind(eval_gauss_beam<Receiver>, receiver_beam, _1));
    auto b = std::chrono::high_resolution_clock::now();
    evalt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    double max_eval_distance_squared = std::pow(max_eval_distance, 2);
    for (size_t source_index = 0; source_index < data.sources().size(); ++source_index) {
        if (squared_distance(source_beam.last_position(), data.sources()[source_index]) >
            max_eval_distance_squared) {
            continue;
        }
        for (size_t receiver_index = 0; receiver_index < data.receivers().size();
             ++receiver_index) {
            if (squared_distance(receiver_beam.last_position(), data.receivers()[receiver_index]) >
                max_eval_distance_squared) {
                continue;
            }
            double total_traveltime =
                std::real(source_beam_values[source_index].complex_traveltime +
                          receiver_beam_values[receiver_index].complex_traveltime);
            if (total_traveltime + window_length > data.timestep() * data.num_samples()) {
                //                std::cerr << total_traveltime << "\n";
                continue;
            }
            a = std::chrono::high_resolution_clock::now();
            Seismogram seismogram = data.get_seismogram(
                data.sources()[source_index], data.receivers()[receiver_index],
                total_traveltime - window_length / 2, total_traveltime + window_length / 2);
            b = std::chrono::high_resolution_clock::now();
            cutt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
            auto seismogram_freq =
                math::fft_closest_frequency(seismogram.data.begin(), seismogram.data.end(),
                                            receiver_beam.frequency(), data.sampling_frequency());
            auto c = std::chrono::high_resolution_clock::now();
            fftt += std::chrono::duration_cast<std::chrono::nanoseconds>(c - b).count();
            stacking_result += source_beam_values[source_index].gb_value *
                               receiver_beam_values[receiver_index].gb_value * seismogram_freq;
        }
    }
    return stacking_result;
}

bool isfinite(const std::complex<double>& c) {
    return std::isfinite(c.real()) & std::isfinite(c.imag());
}

DoubleBeamResult DoubleBeam::algorithm(std::vector<position_t> source_geometry,
                                       std::vector<position_t> target_geometry, SeismoData data,
                                       FractureParameters fracture_info, Meter beam_width,
                                       AngularFrequency beam_frequency, double window_length,
                                       double max_stacking_distance) {
    DoubleBeamResult result(fracture_info.spacings.size(), fracture_info.orientations.size());
    auto ray_code = direct_ray_code(target_geometry[0], source_geometry[0], model);
    std::cout << target_geometry.size() << "Targets\n"
              << source_geometry.size() << "Sources." << std::endl;
    for (const auto& target : target_geometry) {
        int source_beam_index = 0;
        for (const auto& source_beam_center : source_geometry) {
            std::cout << "Beam " << source_beam_index++ << std::endl;
            slowness_t slowness = twopoint.trace(target, source_beam_center);
            // TODO add overload so declaring initial state is not required for ray tracing
            auto initial_state = make_state(target, slowness);
            auto a = std::chrono::high_resolution_clock::now();
            auto source_beam =
                tracer.trace_beam(initial_state, beam_width, beam_frequency, ray_code, {});
            auto b = std::chrono::high_resolution_clock::now();
            beamt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
            if (source_beam.status == Status::OutOfBounds) {
                continue;
            }
            // Since we traced the beam from the target upwards to the surface, we will have to flip
            // the direction of the slowness to be able to treat it as the incoming direction of the
            // beam at the fractures and then scatter.
            std::get<0>(slowness) *= -1;
            std::get<1>(slowness) *= -1;
            std::get<2>(slowness) *= -1;
            for (auto spacing_index = 0U; spacing_index < fracture_info.spacings.size();
                 ++spacing_index) {
                for (auto orientations_index = 0U;
                     orientations_index < fracture_info.orientations.size(); ++orientations_index) {
//                    fmt::print("Spacing {}, orientation {}\n", spacing_index, orientations_index);
                    auto [phi_hat_x, phi_hat_y] = fracture_info.orientations[orientations_index];
                    auto [px, py] = scattered_slowness(
                        std::get<0>(slowness), std::get<1>(slowness), phi_hat_x, phi_hat_y,
                        fracture_info.spacings[spacing_index], angular_to_hertz(beam_frequency));
                    // trace receiver beam in scattered direction
                    // -pz to reflect beam upwards from target
                    slowness_t new_slowness = math::scale_vector(
                        {px, py, -std::get<2>(slowness)}, std::apply(math::length, slowness));
                    initial_state = make_state(target, new_slowness);
                    // reuse ray code since beam should pass through the same layers
                    a = std::chrono::high_resolution_clock::now();
                    auto receiver_beam =
                        tracer.trace_beam(initial_state, beam_width, beam_frequency, ray_code, {});
                    b = std::chrono::high_resolution_clock::now();
                    beamt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
                    if (receiver_beam.status == Status::OutOfBounds) {
                        // beam didn't reach surface, skip
                        std::cout << "Receiver beam left model" << std::endl;
                        continue;
                    }
                    // iteration over sources and receivers
                    auto tmp = stack(source_beam.value(), receiver_beam.value(), data,
                                     window_length, max_stacking_distance);
                    if (not isfinite(tmp)) {
                        std::cerr << "(" << spacing_index << ", " << orientations_index
                                  << ") = " << tmp << std::endl;
                    }
                    result.data(spacing_index, orientations_index) += tmp;
                }
            }
        }
    }
    std::cout << "Beams: " << beamt * 1E-9 << " s\nFFT: " << fftt * 1E-9
              << " s\nBeam eval: " << evalt * 1E-9 << " s\ncut : " << cutt * 1E-9 << "s\n";
    std::cout << "amplitude: " << amplt * 1E-9 << " s\nexp: " << expt * 1E-9
              << " s\nunit vec: " << unit_vect * 1E-9 << " s\nrestt: " << restt * 1E-9;
    return result;
}
