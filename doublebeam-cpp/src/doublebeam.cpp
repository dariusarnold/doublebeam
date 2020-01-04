#include <complex>

#include <boost/range/adaptor/indexed.hpp>
#include <fmt/color.h>
#include <fmt/format.h>

#include "doublebeam.hpp"
#include "eigen_helpers.hpp"
#include "printing.hpp"


/**
 * Calculate horizontal slowness for wave scattered from fractures.
 * @param slow Slowness vector of incoming ray
 * @param phi_hat Unit vector normal to fracture plane. Since only vertical fracture planes are
 * considered, only x and y component of the vector is required.
 * @param fracture_spacing Distance between fracture planes.
 * @param frequency Frequency of wave.
 * @return X, Y component of slowness.
 */
std::tuple<InverseVelocity, InverseVelocity> scattered_slowness(const Slowness& slow,
                                                                const math::Vector2& phi_hat,
                                                                Meter fracture_spacing,
                                                                Frequency frequency) {
    // pass 0 as pz and phi_hat_z since formula only transforms horizontal slownesses and is only
    // valid for vertical fracture planes.
    auto sig = math::sign(math::dot(slow.px.get(), slow.py.get(), 0, phi_hat.x, phi_hat.y, 0));
    InverseVelocity px_new(slow.px.get() -
                           sig * phi_hat.x / (fracture_spacing.get() * frequency.get()));
    InverseVelocity py_new(slow.py.get() -
                           sig * phi_hat.y / (fracture_spacing.get() * frequency.get()));
    return {px_new, py_new};
}


/**
 * Calculate new slowness by scattering and the normalizing length of resulting slowness vector
 * to incoming slowness vector.
 * @param slowness Slowness vector of incoming ray
 * @param fracture_normal Unit vector normal to fracture plane. Since only vertical fracture planes
 * are considered, only x and y component of the vector is required.
 * @param fracture_spacing Distance between fracture planes.
 * @param frequency Frequency of wave.
 * @return Scattered slowness vector
 */
Slowness calculate_new_slowness(const Slowness& slowness, const math::Vector2& fracture_normal,
                                Meter fracture_spacing, Frequency frequency) {
    const auto [px, py] =
        scattered_slowness(slowness, fracture_normal, fracture_spacing, frequency);
    auto p = math::length(slowness);
    // scale px and py so that the whole slowness vector has length 1/c while keeping pz constant
    // p^2 = a*px^2 + a*py^2 + pz^2     where a is the scaling factor
    // a = (p^2 - pz^2) / (px^2 + py^2)
    const double length_factor = std::sqrt((std::pow(p, 2) - std::pow(slowness.pz.get(), 2)) /
                                           (px.get() * px.get() + py.get() * py.get()));
    // -pz to reflect beam upwards from target
    return {px * length_factor, py * length_factor, -slowness.pz};
}


FractureParameters::FractureParameters(math::Vector2 phi_hat_, int num_fracture_orientations,
                                       Meter spacing_min, Meter spacing_max,
                                       int num_fracture_spacings) :
        phi_hat(phi_hat_),
        orientations(math::generate_vector_arc(num_fracture_orientations, phi_hat)),
        spacings(math::linspace(spacing_min, spacing_max, num_fracture_spacings)) {}

std::vector<WaveType> direct_ray_code(Position source, Position receiver,
                                      const VelocityModel& model) {
    auto n = model.number_of_interfaces_between(source.z, receiver.z);
    return std::vector<WaveType>(n, WaveType::Transmitted);
}

DoubleBeam::DoubleBeam(const VelocityModel& model_) :
        model(model_), twopoint(model_), tracer(model_) {}

#define USEDEBUG false
#if USEDEBUG
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

/**
 * Calculate beam amplitude at s.
 * @param beam
 * @param s
 * @return Normalized amplitde of Gauss beam.
 */
std::complex<double> gb_amplitude(const Beam& beam, Arclength s) {
    std::complex<double> det_Q_s = beam.get_Q(s).determinant();
    std::complex<double> det_Q_s0 = beam.get_Q(Arclength(0_meter)).determinant();
    Velocity vs = beam.velocity(s);
    Velocity vs0 = beam.velocity(Arclength(0_meter));
    return std::sqrt((vs.get() * det_Q_s0) / (vs0.get() * det_Q_s));
}

std::complex<double> gb_exp(const Beam& beam, double q1, double q2, Arclength s,
                            std::complex<double>& complex_traveltime) {
    using namespace std::complex_literals;
    Eigen::Vector2d q{q1, q2};
    auto Q = beam.get_Q(s);
    auto M_s = beam.get_P(s) * Q.inverse();
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


Position operator-(const Position& a, const Position& b) {
    return Position(a.x - b.x, a.y - b.y, a.z - b.z);
}

position_t make_tuple(const Position& position) {
    return {position.x.get(), position.y.get(), position.z.get()};
}

BeamEvalResult eval_gauss_beam(const Beam& beam, const Position& position) {
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
    auto [q1, q2, q3] =
        math::dot(transformation_matrix, make_tuple(position - beam.last_position()));
    const Meter total_arclength = beam.last_arclength().length + Meter(q3);
    if (total_arclength < 0_meter) {
        throw std::runtime_error("Negative arclength");
    }
    auto amp = gb_amplitude(beam, Arclength(total_arclength));
    std::complex<double> complex_traveltime;
    auto exp = gb_exp(beam, q1, q2, Arclength(total_arclength), complex_traveltime);
    return {std::conj(amp * exp), complex_traveltime};
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

std::unordered_map<Seismogram<const double>, std::complex<double>,
                   boost::hash<Seismogram<const double>>>
    fft_cache;

std::complex<double> stack(const Beam& source_beam, const Beam& receiver_beam,
                           const SeismoData& data, Second window_length, Meter max_eval_distance) {
    std::complex<double> stacking_result(0, 0);
    std::vector<BeamEvalResult> source_beam_values, receiver_beam_values;
    // get all sources/receivers within max eval distance around surface point
    auto a = std::chrono::high_resolution_clock::now();
    auto sources_in_range = data.get_sources(source_beam.last_position(), max_eval_distance);
    auto receivers_in_range = data.get_receivers(receiver_beam.last_position(), max_eval_distance);
    source_beam_values.reserve(std::size(sources_in_range));
    receiver_beam_values.reserve(std::size(receivers_in_range));
    // pre compute gauss beam for sources/receivers to avoid multiple evaluation
    std::transform(sources_in_range.begin(), sources_in_range.end(),
                   std::back_inserter(source_beam_values),
                   [&](const Source& source) { return eval_gauss_beam(source_beam, source); });
    std::transform(receivers_in_range.begin(), receivers_in_range.end(),
                   std::back_inserter(receiver_beam_values), [&](const Receiver& receiver) {
                       return eval_gauss_beam(receiver_beam, receiver);
                   });
    auto b = std::chrono::high_resolution_clock::now();
    evalt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
    namespace ba = boost::adaptors;
    for (const auto& source : sources_in_range | ba::indexed()) {
        for (const auto& receiver : receivers_in_range | ba::indexed()) {
            Second total_traveltime(
                std::real(source_beam_values[source.index()].complex_traveltime +
                          receiver_beam_values[receiver.index()].complex_traveltime));
            if (total_traveltime - window_length / 2 > data.time_length()) {
                // evaluation position to far away from beam surface point
                continue;
            }
            a = std::chrono::high_resolution_clock::now();
            Seismogram seismogram = data.get_seismogram(source.value(), receiver.value(),
                                                        total_traveltime - window_length / 2,
                                                        total_traveltime + window_length / 2);
            b = std::chrono::high_resolution_clock::now();
            cutt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
            auto it = fft_cache.find(seismogram);
            if (it == fft_cache.end()) {
                auto seismogram_freq = math::fft_closest_frequency(
                    seismogram.data, receiver_beam.frequency(), data.sampling_frequency());
                fft_cache[seismogram] = seismogram_freq;
            }
            auto c = std::chrono::high_resolution_clock::now();
            fftt += std::chrono::duration_cast<std::chrono::nanoseconds>(c - b).count();
            stacking_result += source_beam_values[source.index()].gb_value *
                               receiver_beam_values[receiver.index()].gb_value *
                               fft_cache[seismogram];
        }
    }
    return stacking_result;
}

bool isfinite(const std::complex<double>& c) {
    return std::isfinite(c.real()) and std::isfinite(c.imag());
}


DoubleBeamResult DoubleBeam::algorithm(const std::vector<Position>& source_geometry,
                                       Position target, const SeismoData& data,
                                       FractureParameters fracture_info, Meter beam_width,
                                       AngularFrequency beam_frequency, Second window_length,
                                       Meter max_stacking_distance) {
    DoubleBeamResult result(fracture_info.spacings.size(), fracture_info.orientations.size());
    auto ray_code = direct_ray_code(target, source_geometry[0], model);
    int source_beam_index = 1;
    for (const auto& source_beam_center : source_geometry) {
        fmt::print(fmt::emphasis::bold, "{}/{} source beam centers\n", source_beam_index++,
                   source_geometry.size());
        Slowness slowness = twopoint.trace(target, source_beam_center);
        // TODO add overload so declaring initial state is not required for ray tracing
        auto initial_state = make_state(target, slowness);
        auto a = std::chrono::high_resolution_clock::now();
        auto source_beam =
            tracer.trace_beam(initial_state, beam_width, beam_frequency, ray_code, {});
        auto b = std::chrono::high_resolution_clock::now();
        beamt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
        if (source_beam.status == Status::OutOfBounds) {
            throw std::logic_error("Source beam left model.\n");
        }
        // Since we traced the beam from the target upwards to the surface, we will have to flip
        // the direction of the slowness to be able to treat it as the incoming direction of the
        // beam at the fractures and then scatter.
        slowness.flip_direction();
        int number_of_rec_beams_that_left_model = 0;
        namespace ba = boost::adaptors;
        for (const auto& fracture_spacing : fracture_info.spacings | ba::indexed(0)) {
            for (const auto& fracture_orientation : fracture_info.orientations | ba::indexed(0)) {
                //                    fmt::print("Spacing {}, orientation {}\n", spacing_index,
                //                    orientations_index);
                // trace receiver beam in scattered direction
                Slowness new_slowness = calculate_new_slowness(
                    slowness, fracture_orientation.value(), fracture_spacing.value(),
                    angular_to_hertz(beam_frequency));
                initial_state = make_state(target, new_slowness);
                // reuse ray code since beam should pass through the same layers
                a = std::chrono::high_resolution_clock::now();
                auto receiver_beam =
                    tracer.trace_beam(initial_state, beam_width, beam_frequency, ray_code);
                b = std::chrono::high_resolution_clock::now();
                beamt += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
                if (receiver_beam.status == Status::OutOfBounds) {
                    // beam didn't reach surface, skip
                    number_of_rec_beams_that_left_model++;
                    continue;
                }
                // iteration over sources and receivers
                auto tmp = stack(source_beam.value(), receiver_beam.value(), data, window_length,
                                 max_stacking_distance);
                if (not isfinite(tmp)) {
                    throw(std::runtime_error(impl::Formatter()
                                             << "(" << fracture_spacing.index() << ", "
                                             << fracture_orientation.index() << ") = " << tmp));
                }
                result.data(fracture_spacing.index(), fracture_orientation.index()) += tmp;
            }
        }
        fmt::print("{}/{} receiver beams left the model.\n", number_of_rec_beams_that_left_model,
                   fracture_info.spacings.size() * fracture_info.orientations.size());
    }
    fmt::print(
        "Beams: {} s\nFFT: {} s\nBeam eval: {} s\ncuting seismograms: {} s\nGB amplitude: {} "
        "s\nGB exp: {} s\nUnit vectors: {} s\nRest: {} s\n",
        beamt * 1E-9, fftt * 1E-9, evalt * 1E-9, cutt * 1E-9, amplt * 1E-9, expt * 1E-9,
        unit_vect * 1E-9, restt * 1E-9);
    return result;
}
