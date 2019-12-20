#include "raytracing.hpp"

#include <Eigen/Dense>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <cstdint>
#include <unsupported/Eigen/CXX11/Tensor>

#include "eigen_helpers.hpp"
#include "model.hpp"
#include "ray.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"


RayState init_state(Meter x, Meter y, Meter z, const VelocityModel& model, Radian theta, Radian phi,
                    TravelTime T) {
    if (not model.in_model(x, y, z)) {
        throw std::domain_error(impl::Formatter() << "Point (" << x << ", " << y << ", " << z
                                                  << ") not in model " << model);
    }
    Velocity velocity(model.eval_at(x, y, z).value());
    auto slowness = seismo::slowness_3D(theta, phi, velocity);
    return {Position{x, y, z}, slowness, TravelTime{T}, Arclength{0_meter}};
}

RayState init_state(Position position, const VelocityModel& model, Radian theta, Radian phi,
                    TravelTime T) {
    return init_state(position.x, position.y, position.z, model, theta, phi, T);
}

RayState make_state(Meter x, Meter y, Meter z, InverseVelocity px, InverseVelocity py,
                    InverseVelocity pz, TravelTime T) {
    return {Position{x, y, z}, Slowness{px, py, pz}, T, Arclength{0._meter}};
}

RayState make_state(Position position, Slowness slowness, TravelTime T) {
    return {position, slowness, T, Arclength{0_meter}};
}

RayTracingResult<RaySegment> RayTracer::trace_layer(const RayState& initial_state,
                                                    const Layer& layer) {
    auto [position, slowness, time, arclength] = initial_state;
    const auto c = layer.velocity;
    bool stop_depth_was_reached = false;
    const Meter z_end = [&]() {
        // check if stop depth exists and we are in the layer where we have to stop
        if (stop_depth_m and
            math::between(layer.top_depth, stop_depth_m.value(), layer.bot_depth)) {
            stop_depth_was_reached = true;
            return stop_depth_m.value();
        } else {
            return seismo::ray_direction_down(slowness) ? layer.bot_depth : layer.top_depth;
        }
    }();
    auto arclength_in_layer = (z_end.get() - position.z.get()) / (c.get() * slowness.pz.get());
    // s has to start at 0 because this calculation is done starting from the current point of
    // the ray.
    const Meter x_end(position.x.get() + arclength_in_layer * c.get() * slowness.px.get());
    const Meter y_end(position.y.get() + arclength_in_layer * c.get() * slowness.py.get());
    if (not model.in_horizontal_extent(x_end, y_end)) {
        // ray left model to the side and did not reach top/bottom of layer

        return {Status::OutOfBounds, {}};
    }
    const Second t_end(time.time.get() + arclength_in_layer / layer.velocity.get());
    const Meter s_end(arclength.length.get() + arclength_in_layer);
    return {stop_depth_was_reached ? Status::StopDepthReached : Status::Success,
            RaySegment{initial_state,
                       RayState{Position{x_end, y_end, z_end}, Slowness{slowness},
                                TravelTime{t_end}, Arclength{s_end}},
                       Velocity(layer.velocity)}};
}


RayTracer::RayTracer(VelocityModel velocity_model) : model(std::move(velocity_model)) {}

RayTracingResult<Ray> RayTracer::trace_ray(const RayState& initial_state,
                                           const std::vector<WaveType>& ray_code,
                                           std::optional<Meter> stop_depth) {
    if (not model.in_model(initial_state.position.x, initial_state.position.y,
                           initial_state.position.z)) {
        throw std::domain_error(impl::Formatter()
                                << "Point " << initial_state.position << " not in model " << model);
    }
    stop_depth_m = stop_depth;
    int64_t layer_index = model.layer_index(initial_state.position.z).value();
    auto current_layer = model[layer_index];
    auto segment = trace_layer(initial_state, current_layer);
    if (segment.status == Status::OutOfBounds) {
        // early exit when ray left model horizontally
        return {Status::OutOfBounds, {}};
    }
    if (segment.status == Status::StopDepthReached) {
        return {Status::StopDepthReached, Ray{{segment.value()}}};
    }
    if (ray_code.empty()) {
        // only trace one layer for empty ray code
        return {Status::Success, {{Ray{segment.value()}}}};
    }
    Ray ray{{segment.value()}};
    for (auto ray_type : ray_code) {
        // TODO this line and the trace_layer call below violates the law of Demeter, try to
        //  refactor it by improving Ray class
        const auto& last_state = ray.last_state();
        // reflected waves stay in the same layer, so the index doesn't change
        if (ray_type == WaveType::Transmitted) {
            layer_index += seismo::ray_direction_down(last_state.slowness.pz.get()) ? 1 : -1;
            current_layer = model[layer_index];
        }
        if (layer_index < 0 or layer_index >= static_cast<int64_t>(model.num_layers())) {
            throw std::runtime_error(impl::Formatter() << "Ray left model at top or bottom.");
        }
        auto new_slowness = snells_law(last_state, model, ray_type);
        RayState new_initial_state{last_state.position, new_slowness, last_state.travel_time,
                                   last_state.arclength};
        segment = trace_layer(new_initial_state, current_layer);
        if (segment.status == Status::OutOfBounds) {
            return {Status::OutOfBounds, {}};
        }
        ray.add_segment(segment.value());
        if (segment.status == Status::StopDepthReached) {
            return {Status::StopDepthReached, ray};
        }
    }
    return {Status::Success, ray};
}


#define USEDEBUG false
#if USEDEBUG
#define msg(x) std::cerr << #x << ": " << x << std::endl;
#else
#define msg(x)
#endif


class InterfacePropagator {
    using matrix_t = Eigen::Matrix2cd;

public:
    /**
     * Transform matrix P and Q across an interface in a velocity model of horizontal layers.
     * Procedure is described in Seismic Ray Theory - Cerveny 2001.
     * @param P Value of P (shape 2x2) at the interface.
     * @param Q Value of Q (shape 2x2) at the interface).
     * @param wave_type Wave type, valid 'R' for reflected or 'T' for transmitted.
     * @param old_state State before interface crossing (position, slowness, trave time).
     * @param new_state State after interface crossing (position, slowness, trave time).
     * @param layer_index Index of the layer the wave is in before the interface is crossed.
     * @param model Velocity model.
     * @return New values for P, Q.
     */
    std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd>
    transform(const Eigen::Matrix2cd& P, const Eigen::Matrix2cd& Q, WaveType wave_type,
              const RayState& old_state, const RayState& new_state, const VelocityModel& model) {
        // TODO modify interface unit vector (params x2, y2, z2) for more general velocity model.
        //  Here it is assumed the model consists only of horizontal layers.
        msg(layer_index);
        msg(wave_type);
        auto i_S = math::angle(old_state.slowness.px.get(), old_state.slowness.py.get(),
                               old_state.slowness.pz.get(), 0, 0, 1);
        msg(i_S);
        auto i_R = wave_type == WaveType::Transmitted
                       ? math::angle(new_state.slowness.px.get(), new_state.slowness.py.get(),
                                     new_state.slowness.pz.get(), 0, 0, 1)
                       : i_S;
        msg(i_R);
        // epsilon is introduced by eq. 2.4.71, Cerveny2001. This formula is simplified for
        // horizontal interfaces (unit vector (0, 0, 1)).
        auto epsilon = std::copysign(1., old_state.slowness.pz.get());
        msg(epsilon);
        // for a downgoing transmitted ray the velocity above the interface is the before
        // velocity and the velocity below the interface is the after velocity.
        auto [V_top, V_bottom] = model.interface_velocities(old_state.position.z);
        auto V_before = V_top, V_after = V_bottom;
        if (wave_type == WaveType::Reflected) {
            V_after = V_before;
        } else {
            if (not seismo::ray_direction_down(old_state.slowness.pz.get())) {
                std::swap(V_before, V_after);
            }
        }
        msg(V_after);
        msg(V_before);
        // TODO this kappa is only valid for simple velocity model v = v(z) and horizontal
        //  interfaces
        auto kappa = 0.;
        auto cos_kappa = std::cos(kappa), sin_kappa = std::sin(kappa);
        // right equations of (4.4.49) in Cerveny2001
        matrix_t G_orthogonal;
        G_orthogonal << cos_kappa, -sin_kappa, sin_kappa, cos_kappa;
        msg(G_orthogonal);
        auto G_orthogonal_tilde = G_orthogonal;
        // left equations of (4.4.49) in Cerveny2001
        matrix_t G_parallel;
        G_parallel << epsilon * std::cos(i_S), 0, 0, 1;
        msg(G_parallel);
        matrix_t G_parallel_tilde;
        G_parallel_tilde << (wave_type == WaveType::Transmitted ? 1 : -1) * epsilon * std::cos(i_R),
            0, 0, 1;
        msg(G_parallel_tilde);
        // equation (4.4.48) from Cerveny2001
        auto G = G_parallel * G_orthogonal;
        msg(G);
        matrix_t G_tilde = G_parallel_tilde * G_orthogonal_tilde;
        msg(G_tilde);
        // Evaluate this since it is used two times and would be reevaluated otherwise
        matrix_t G_inverted = G.inverse();
        msg(G_inverted);
        // TODO simplify by adapting for constant velocity layers
        auto old_gradient = 0;
        msg(old_gradient);
        auto new_gradient = 0;
        msg(new_gradient);
        // eq. (4.4.53) from Cerveny2001
        auto E = E_(V_before.get(), i_S, epsilon, old_gradient);
        msg(E);
        auto E_tilde = E_tilde_(wave_type, V_after.get(), i_R, epsilon, new_gradient);
        msg(E_tilde);
        auto u = u_(wave_type, V_before.get(), V_after.get(), i_S, i_R, epsilon);
        msg(u);
        auto D = D_();
        // eq. (4.4.67) Cerveny2001
        msg(D);
        matrix_t P_tilde =
            G_tilde.inverse() * ((G * P) + (E - E_tilde - u * D) * (G_inverted.transpose() * Q));
        // eq. (4.4.64) from Cerveny2001
        msg(P_tilde);
        matrix_t Q_tilde = G_tilde.transpose() * G_inverted.transpose() * Q;
        msg(Q_tilde);
        return {P_tilde, Q_tilde};
    }

private:
    /**
     * Eq. 4.4.53 from Cerveny2001
     * @param V
     * @param i_S
     * @param epsilon
     * @param old_gradient
     * @return
     */
    matrix_t E_(double V, double i_S, double epsilon, double old_gradient) const {
        // TODO modify this to work with a more general velocity model
        // dV_dzi means the derivative of the velocity after the z_i coordinate for V=V(z)
        auto dV_dz1 = 0.;
        auto dV_dz2 = 0.;
        auto dV_dz3 = old_gradient;
        auto E11 = -std::sin(i_S) / (V * V) *
                   ((1 + std::pow(std::cos(i_S), 2)) * dV_dz1 -
                    epsilon * std::cos(i_S) * std::sin(i_S) * dV_dz3);
        auto E12 = -std::sin(i_S) / (V * V) * dV_dz2;
        auto E22 = 0.;
        return (matrix_t() << E11, E12, E12, E22).finished();
    }

    /**
     * Eq. 4.4.54 from Cerveny2001
     * @param wave_type
     * @param V_tilde
     * @param i_R
     * @param epsilon
     * @param new_gradient
     * @return
     */
    matrix_t E_tilde_(WaveType wave_type, double V_tilde, double i_R, double epsilon,
                      double new_gradient) const {
        auto dV_tilde_dz1 = 0.;
        auto dV_tilde_dz2 = 0.;
        auto dV_tilde_dz3 = new_gradient;
        auto minus_plus = wave_type == WaveType::Reflected ? -1. : 1.;
        auto E11 = -std::sin(i_R) / (V_tilde * V_tilde) *
                   ((1 + std::pow(cos(i_R), 2)) * dV_tilde_dz1 +
                    minus_plus * epsilon * std::cos(i_R) * std::sin(i_R) * dV_tilde_dz3);
        auto E12 = -std::sin(i_R) / (V_tilde * V_tilde) * dV_tilde_dz2;
        auto E22 = 0.;
        return (matrix_t() << E11, E12, E12, E22).finished();
    }

    /**
     * Eq. 4.4.51 from Cerveny2001
     * @param wave_type String specifying wave type, valid values are "T" for
        transmitted and "R" for reflected.
     * @param V Velocity before the interface, in m/s.
     * @param V_tilde Velocity after the interface, in m/s.
     * @param i_S Acute angle of incidence, 0 <= i_S <= pi/2.
     * @param i_R Acute angle of reflection/transmission
     * @param epsilon sign(p * n)
     */
    static double u_(WaveType wave_type, double V, double V_tilde, double i_S, double i_R,
                     double epsilon) {
        // TODO make this a function
        auto minusplus = wave_type == WaveType::Reflected ? -1 : 1;
        return epsilon * (std::cos(i_S) / V + minusplus * std::cos(i_R) / V_tilde);
    }

    /**
     * Eq. 4.4.15 from Cerveny2001
     * For the currently implemented velocity layer with horizontal interfaces only, this function
     * is zero everywhere since it contains the second derivative of the interface function Sigma
     * in the numerator. Sigma = Sigma(z3) for horizontal interfaces.
     */
    static matrix_t D_() {
        return matrix_t::Zero();
    }
};


/**
 * Calculate distance between point a and b.
 */
Meter distance(const Position& a, const Position& b) {
    auto [ax, ay, az] = a;
    auto [bx, by, bz] = b;
    return Meter(std::sqrt(std::pow(ax.get() - bx.get(), 2) + std::pow(ay.get() - by.get(), 2) +
                           std::pow(az.get() - bz.get(), 2)));
}

bool not_at_last_ray_segment(size_t segment_index, size_t ray_size) {
    return segment_index < ray_size - 1;
}

RayTracingResult<Beam> RayTracer::trace_beam(const RayState& initial_state, Meter beam_width,
                                             AngularFrequency beam_frequency,
                                             const std::vector<WaveType>& ray_code,
                                             std::optional<Meter> stop_depth) {
    // first trace ray kinematically
    auto ray = trace_ray(initial_state, ray_code, stop_depth);
    if (not ray.result) {
        // ray tracing failed
        return {ray.status, {}};
    }
    InterfacePropagator ip;
    auto [position, slowness, traveltime, arclength] = initial_state;
    auto v0 = model.eval_at(position).value();
    // initial values for P, Q
    Eigen::Matrix2cd P;
    P << 1j / v0.get(), 0, 0, 1j / v0.get();
    Eigen::Matrix2cd Q;
    Q << beam_frequency.get() * beam_width.get() * beam_width.get() / v0.get(), 0, 0,
        beam_frequency.get() * beam_width.get() * beam_width.get() / v0.get();
    Beam beam(beam_width, beam_frequency);
    std::ptrdiff_t segment_index = 0;
    for (const auto& segment : ray.value()) {
        beam.add_segment(P, Q, segment);
        if (not_at_last_ray_segment(segment_index, ray.value().size())) {
            // if we are not at the last segment of the ray, transform dynamic ray tracing across
            // interface (calculate new P0, Q0)
            Eigen::Matrix2cd Q_at_end_of_ray_segment =
                Q + segment.layer_velocity().get() * segment.length().get() * P;
            auto wave_type = ray_code[segment_index];
            auto new_initial_state = ray.value()[segment_index + 1].begin();
            std::tie(P, Q) = ip.transform(P, Q_at_end_of_ray_segment, wave_type, segment.end(),
                                          new_initial_state, model);
        }
        ++segment_index;
    }
    return {Status::Success, beam};
}

RayTracingResult<Ray> RayTracer::trace_ray(const RayState& initial_state,
                                           const std::string& ray_code,
                                           std::optional<Meter> stop_depth) {
    return trace_ray(initial_state, seismo::make_ray_code(ray_code), stop_depth);
}

RayTracingResult<Beam> RayTracer::trace_beam(const RayState& initial_state, Meter beam_width,
                                             AngularFrequency beam_frequency,
                                             const std::string& ray_code,
                                             std::optional<Meter> stop_depth) {
    return trace_beam(initial_state, beam_width, beam_frequency, seismo::make_ray_code(ray_code),
                      stop_depth);
}
