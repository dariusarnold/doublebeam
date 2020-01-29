#include <complex>

#include <boost/range/adaptor/indexed.hpp>
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

struct UnitVectors {
    using vector_t = Eigen::RowVector3d;
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
    const UnitVectors::vector_t z_axis = (UnitVectors::vector_t() << 0, 0, 1).finished();
    Eigen::Vector3d e2 = z_axis.cross(Eigen::Vector3d(px.get(), py.get(), 0));
    // find other unit vector
    const UnitVectors::vector_t p(px.get(), py.get(), pz.get());
    UnitVectors::vector_t e1 = -1 * p.cross(e2);
    e1.normalize();
    e2.normalize();
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

std::complex<double> complex_traveltime(const Beam& beam, double q1, double q2, Arclength s) {
    Eigen::Vector2d q{q1, q2};
    Eigen::Matrix2cd Q = beam.get_Q(s);
    return beam.traveltime().get() + 0.5 * (q.transpose() * beam.get_P(s) * Q.inverse() * q)[0];
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


struct BeamEvalResult {
    std::complex<double> gb_value;
    std::complex<double> amplitude;
    std::complex<double> complex_traveltime;
};

Position operator-(const Position& a, const Position& b) {
    return Position(a.x - b.x, a.y - b.y, a.z - b.z);
}

Eigen::Vector3d make_vector(const Position& position) {
    return {position.x.get(), position.y.get(), position.z.get()};
}

BeamEvalResult eval_gauss_beam(const Beam& beam, const Position& position, AngularFrequency omega) {
    using namespace std::complex_literals;
    auto [e1, e2] = get_ray_centred_unit_vectors(beam);
    UnitVectors::vector_t e3 = e1.cross(e2);
    e3.normalize();
    // We want to convert from general Cartesian system to ray centred, use right eq. 4.1.31 from
    // Cerveny2001.
    Eigen::Matrix3d transformation_matrix;
    transformation_matrix << e1, e2, e3;
    const Eigen::Vector3d q = transformation_matrix * make_vector(position - beam.last_position());
    const Meter total_arclength = beam.last_arclength().length + Meter(q[2]);
    if (total_arclength < 0_meter) {
        throw std::runtime_error("Negative arclength");
    }
    auto amp = gb_amplitude(beam, Arclength(total_arclength));
    std::complex<double> traveltime =
        complex_traveltime(beam, q[0], q[1], Arclength(total_arclength));
    return {std::conj(amp * std::exp(1i * omega.get() * traveltime)), amp, traveltime};
}


DoubleBeamResult::DoubleBeamResult(size_t num_of_fracture_spacings,
                                   size_t num_of_fracture_orientations) :
        data(num_of_fracture_spacings, num_of_fracture_orientations) {
    data.setZero();
}

std::complex<double> stack(const Beam& source_beam, const Beam& receiver_beam,
                           const SeismoData& data, Second window_length, Meter max_eval_distance,
                           AngularFrequency source_frequency) {
    std::complex<double> stacking_result(0, 0);
    std::vector<BeamEvalResult> source_beam_values, receiver_beam_values;
    // get all sources/receivers within max eval distance around surface point
    auto sources_in_range = data.get_sources(source_beam.last_position(), max_eval_distance);
    auto receivers_in_range = data.get_receivers(receiver_beam.last_position(), max_eval_distance);
    source_beam_values.reserve(std::size(sources_in_range));
    receiver_beam_values.reserve(std::size(receivers_in_range));
    // pre compute gauss beam for sources/receivers to avoid multiple evaluation
    std::transform(sources_in_range.begin(), sources_in_range.end(),
                   std::back_inserter(source_beam_values), [&](const Source& source) {
                       return eval_gauss_beam(source_beam, source, source_frequency);
                   });
    std::transform(receivers_in_range.begin(), receivers_in_range.end(),
                   std::back_inserter(receiver_beam_values), [&](const Receiver& receiver) {
                       return eval_gauss_beam(receiver_beam, receiver, source_frequency);
                   });
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
            Seismogram seismogram_x = data.get_seismogram<Component::X>(
                source.value(), receiver.value(), total_traveltime - window_length / 2,
                total_traveltime + window_length / 2);
            Seismogram seismogram_y = data.get_seismogram<Component::Y>(
                source.value(), receiver.value(), total_traveltime - window_length / 2,
                total_traveltime + window_length / 2);
            Seismogram seismogram_z = data.get_seismogram<Component::Z>(
                source.value(), receiver.value(), total_traveltime - window_length / 2,
                total_traveltime + window_length / 2);
            auto seismogram_freq_x =
                math::fft(seismogram_x.data, receiver_beam.frequency(), data.sampling_frequency());
            auto seismogram_freq_y =
                math::fft(seismogram_y.data, receiver_beam.frequency(), data.sampling_frequency());
            auto seismogram_freq_z =
                math::fft(seismogram_z.data, receiver_beam.frequency(), data.sampling_frequency());
            stacking_result += source_beam_values[source.index()].gb_value *
                               receiver_beam_values[receiver.index()].gb_value *
                               (seismogram_freq_x + seismogram_freq_y + seismogram_freq_z) / 3.;
        }
    }
    std::complex<double> value{0, 0};
    for (const auto& receiver : receivers_in_range | ba::indexed()) {
        value += std::pow(std::abs(receiver_beam_values[receiver.index()].amplitude), 2);
    }
    if (receivers_in_range.size() == 0) {
        return stacking_result;
    }
    return stacking_result / value;
}

DoubleBeamResult DoubleBeam::algorithm(const std::vector<Position>& source_geometry,
                                       Position target, const SeismoData& data,
                                       const FractureParameters& fracture_info, Meter beam_width,
                                       AngularFrequency beam_frequency, Second window_length,
                                       Meter max_stacking_distance,
                                       AngularFrequency source_frequency) {
    DoubleBeamResult result(fracture_info.spacings.size(), fracture_info.orientations.size());
    auto ray_code = direct_ray_code(target, source_geometry[0], model);
    int source_beam_index = 0;
    Eigen::ArrayXXcd temp(result.data);
    #pragma omp declare reduction(+:Eigen::ArrayXXcd:omp_out=omp_out+omp_in) \
    initializer(omp_priv=Eigen::ArrayXXcd::Zero(omp_orig.rows(), omp_orig.cols()))
    #pragma omp parallel for reduction(+:temp) schedule(dynamic) default(none)\
    shared(source_geometry, data, fracture_info, source_beam_index, std::cout)\
    firstprivate(beam_width, target, max_stacking_distance, window_length, beam_frequency, ray_code, source_frequency)
    for (auto sbc = source_geometry.begin(); sbc != source_geometry.end(); ++sbc) {
    #pragma omp critical
        {
            fmt::print("\r{}/{} source beam centers", ++source_beam_index, source_geometry.size());
            std::cout.flush();
        }
        temp +=
            calc_sigma_for_sbc(*sbc, target, fracture_info, data, beam_width, beam_frequency,
                               ray_code, window_length, max_stacking_distance, source_frequency);
    }
    result.data = temp;
    // Add newline after the loop progress output
    fmt::print("\n");
    using namespace std::complex_literals;
    result.data *= 2i * source_frequency.get();
    return result;
}

Eigen::ArrayXXcd
DoubleBeam::calc_sigma_for_sbc(const Position& source_beam_center, const Position& target,
                               const FractureParameters& fracture_info, const SeismoData& data,
                               Meter beam_width, AngularFrequency beam_frequency,
                               const std::vector<WaveType>& ray_code, Second window_length,
                               Meter max_stacking_distance, AngularFrequency source_frequency) {
    Eigen::ArrayXXcd result =
        Eigen::ArrayXXcd::Zero(fracture_info.spacings.size(), fracture_info.orientations.size());
    Slowness slowness = twopoint.trace(target, source_beam_center);
    auto source_beam = tracer.trace_beam(target, slowness, beam_width, beam_frequency, ray_code);
    if (source_beam.status == Status::OutOfBounds) {
        throw std::logic_error("Source beam left model.\n");
    }
    // Since we traced the beam from the target upwards to the surface, we will have to flip
    // the direction of the slowness to be able to treat it as the incoming direction of the
    // beam at the fractures and then scatter.
    slowness.flip_direction();
    namespace ba = boost::adaptors;
    for (const auto& fracture_spacing : fracture_info.spacings | ba::indexed()) {
        for (const auto& fracture_orientation : fracture_info.orientations | ba::indexed()) {
            // trace receiver beam in scattered direction
            Slowness new_slowness =
                calculate_new_slowness(slowness, fracture_orientation.value(),
                                       fracture_spacing.value(), angular_to_hertz(beam_frequency));
            // reuse ray code since beam should pass through the same layers
            auto receiver_beam =
                tracer.trace_beam(target, new_slowness, beam_width, beam_frequency, ray_code);
            if (receiver_beam.status == Status::OutOfBounds) {
                // beam didn't reach surface, skip
                continue;
            }
            result(fracture_spacing.index(), fracture_orientation.index()) +=
                stack(source_beam.value(), receiver_beam.value(), data, window_length,
                      max_stacking_distance, source_frequency);
        }
    }
    result *= source_beam.value().last_slowness().pz.get();
    return result;
}
