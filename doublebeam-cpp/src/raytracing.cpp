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
    const auto& slowness = initial_state.slowness;
    const auto& position = initial_state.position;
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
    const Second t_end(initial_state.travel_time.time.get() +
                       arclength_in_layer / layer.velocity.get());
    const Meter s_end(initial_state.arclength.length.get() + arclength_in_layer);
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
#define msg(x) std::cout << #x << ": " << x << std::endl;
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
     * @return New values for P, Q.
     */
    std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd>
    transform(const Eigen::Matrix2cd& P, const Eigen::Matrix2cd& Q, WaveType wave_type,
              const RayState& old_state, const RayState& new_state) {
        // i_S is the acute angle of incidence, 0 <= i_s <= pi/2
        const auto i_S = math::angle(old_state.slowness.px.get(), old_state.slowness.py.get(),
                                     old_state.slowness.pz.get(), 0, 0, 1);
        msg(degrees(i_S));
        // i_r is the acute angle of reflection/transmission
        const auto i_R = wave_type == WaveType::Transmitted
                             ? math::angle(new_state.slowness.px.get(), new_state.slowness.py.get(),
                                           new_state.slowness.pz.get(), 0, 0, 1)
                             : i_S;
        msg(degrees(i_R));
        // epsilon is the orientation index introduced by eq. 2.4.71, Cerveny2001. This formula
        // is simplified for horizontal interfaces (unit vector (0, 0, 1)). Base formula is
        // eps = sign(dot(p, n)) where p is the slowness vector, n is the normal of the interface
        // and sign gives the sign of its argument.
        auto epsilon = std::copysign(1., old_state.slowness.pz.get());
        msg(epsilon);
        // kappa is the angle between e_2 and i_2, 0 <= kappa <= 2pi.
        // i2 is the basis vector of the local Cartesian coordinate system with origin at Q.
        // i3 coincides with the unit vector n normal to the interface, i1 and i2 can be chosen
        // arbitrarily in the plane. In our case i2 will always be horizontal due to the plane
        // interfaces.
        // e2 is a basis vector of the ray centred coordinate system of the incident wave at Q.
        // For the special case of a planar ray e2 will always be horizontal (see function
        // get_ray_centred_unit_vectors). Since we can choose i_1, i_2 so that i is an orthogonal
        // right handed coordinate system, chose i_2 to be parallel to e_2.
        const auto kappa = 0.;
        const auto cos_kappa = std::cos(kappa), sin_kappa = std::sin(kappa);
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
        // simplified version for homogeneous layer with constant velocity and planar horizontal
        // interfaces Cerveny2001 (4.8.10).
        matrix_t P_tilde = G_tilde.inverse() * G * P;
        // eq. (4.4.64) from Cerveny2001
        msg(P_tilde);
        matrix_t Q_tilde = G_tilde.transpose() * G.inverse().transpose() * Q;
        msg(Q_tilde);
        return {P_tilde, Q_tilde};
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
    Q << beam_frequency.get() * std::pow(beam_width.get(), 2) / v0.get(), 0, 0,
        beam_frequency.get() * std::pow(beam_width.get(), 2) / v0.get();
    Beam beam(beam_width, beam_frequency);
    std::ptrdiff_t segment_index = 0;
    for (const auto& segment : ray.value()) {
        beam.add_segment(P, Q, segment);
        if (not_at_last_ray_segment(segment_index, ray.value().size())) {
            // if we are not at the last segment of the ray, transform dynamic ray tracing
            // across interface (calculate new P0, Q0)
            Eigen::Matrix2cd Q_at_end_of_ray_segment =
                Q + segment.layer_velocity().get() * segment.length().get() * P;
            auto wave_type = ray_code[segment_index];
            auto new_initial_state = ray.value()[segment_index + 1].begin();
            std::tie(P, Q) = ip.transform(P, Q_at_end_of_ray_segment, wave_type, segment.end(),
                                          new_initial_state);
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
