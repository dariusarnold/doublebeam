#include "raytracing.hpp"

#include <Eigen/Dense>
#include <boost/numeric/odeint/stepper/generation/generation_dense_output_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include "eigen_helpers.hpp"
#include "model.hpp"
#include "ray.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"


std::string point_to_str(double x, double y, double z) {
    return impl::Formatter(",") << x << y << z;
}

state_type init_state(double x, double y, double z, const VelocityModel& model, double theta,
                      double phi, double T) {
    if (not model.in_model(x, y, z)) {
        throw std::domain_error(impl::Formatter()
                                << "Point " << point_to_str(x, y, z) << " not in model " << model);
    }
    double velocity = model.eval_at(x, y, z).value();
    auto [px, py, pz] = seismo::slowness_3D(theta, phi, velocity);
    return {x, y, z, px, py, pz, T};
}

state_type init_state(position_t position, const VelocityModel& model, double theta, double phi,
                      double T) {
    auto [x, y, z] = position;
    return init_state(x, y, z, model, theta, phi, T);
}

state_type make_state(double x, double y, double z, double px, double py, double pz, double T) {
    return {x, y, z, px, py, pz, T};
}

state_type make_state(position_t position, slowness_t slowness, double T) {
    auto [x, y, z] = position;
    auto [px, py, pz] = slowness;
    return {x, y, z, px, py, pz, T};
}


/**
 * Calculate first derivative of inverse of velocity after depth z analytically.
 * Valid for linear velocity gradient v = v(z) = a * z + b).
 * @param z
 * @param layer
 * @return Derivative d/dz of 1/v(z) = -(az+b)^{-2}*a
 */
inline double dvdz(double z, const Layer& layer) {
    return -layer.gradient /
           ((layer.gradient * z + layer.intercept) * (layer.gradient * z + layer.intercept));
}

/**
 * This function implements the system of ODEs required to compute the ray.
 * The method is not called directly from my code, only by the solver.
 * @param state Current state is input from this.
 * @param dfds Next step is stored here.
 * @param s Current arclength along the ray. The ray tracing system of ODEs does not depend
 * on this parameter.
 */
void ray_tracing_equation(const state_type& state, state_type& dfds, const double /* s */,
                          const VelocityModel& model, const Layer& layer) {
    auto [x, y, z, px, py, pz, T] = state;
    auto v = model.eval_at(x, y, z).value_or(0);
    auto dxds = px * v;
    auto dyds = py * v;
    auto dzds = pz * v;
    auto dpzds = dvdz(z, layer);
    auto dTds = v == 0 ? 0. : 1. / v;
    dfds[0] = dxds;
    dfds[1] = dyds;
    dfds[2] = dzds;
    dfds[3] = 0.;
    dfds[4] = 0.;
    dfds[5] = dpzds;
    dfds[6] = dTds;
}

double layer_height(const Layer& layer) {
    return layer.bot_depth - layer.top_depth;
}

RayTracingResult<RaySegment> RayTracer::trace_layer_gradient(const state_type& initial_state,
                                                             const Layer& layer, double s_start,
                                                             double ds, double max_ds) {
    namespace odeint = boost::numeric::odeint;
    using stepper_t = odeint::runge_kutta_dopri5<state_type>;
    InterfaceCrossed crossing = [&]() {
        if (stop_depth_m) {
            // set stop depth as an interface where integration will stop.
            if (initial_state[Index::Z] > stop_depth_m.value()) {
                // Integration starts at bottom side of layer
                return InterfaceCrossed(stop_depth_m.value(), layer.bot_depth);
            } else {
                return InterfaceCrossed(layer.top_depth, stop_depth_m.value());
            }
        } else {
            // use layer boundaries as integration borders
            return InterfaceCrossed(layer.top_depth, layer.bot_depth);
        }
    }();
    // error values can be lowered to 1e-8 with only minimal loss in precision to improve speed
    auto stepper = odeint::make_dense_output(1.E-10, 1.E-10, max_ds, stepper_t());
    stepper.initialize(initial_state, s_start, ds);
    std::vector<double> arclengths;
    std::vector<state_type> states;
    auto estimate_steps = 2 * layer_height(layer) / ds;
    arclengths.reserve(estimate_steps);
    states.reserve(estimate_steps);
    auto equation = [&](const state_type& state, state_type& dfds, double s) {
        return ray_tracing_equation(state, dfds, s, model, layer);
    };
    // advance stepper until first event occurs
    do {
        if (not model.in_horizontal_extent(stepper.current_state()[Index::X],
                                           stepper.current_state()[Index::Y])) {
            return {Status::OutOfBounds, {}};
        }
        states.emplace_back(stepper.current_state());
        arclengths.emplace_back(stepper.current_time());
        stepper.do_step(equation);
    } while (not crossing(stepper.current_state()));
    // find exact point of zero crossing
    auto crossing_function = crossing.get_zero_crossing_event_function(stepper.current_state());
    auto [state_at_crossing, s_at_crossing] = get_state_at_interface(crossing_function, stepper);
    arclengths.emplace_back(s_at_crossing);
    states.emplace_back(state_at_crossing);
    // workaround for numerical issues where layer boundaries are overshot
    // do this by clamping the value to the interface depth
    states.back()[Index::Z] = crossing.get_closest_layer_depth(states.back());
    bool stop_depth_was_reached = stop_depth_m and states.back()[Index::Z] == stop_depth_m.value();
    return {stop_depth_was_reached ? Status::StopDepthReached : Status::Success,
            RaySegment{states, arclengths}};
}


RayTracingResult<RaySegment> RayTracer::trace_layer_const(const state_type& initial_state,
                                                          const Layer& layer, double s_start,
                                                          double ds) {
    auto [x0, y0, z0, px0, py0, pz0, t0] = initial_state;
    auto c = layer.intercept;
    bool stop_depth_was_reached = false;
    double z_interface = [&]() {
        // check if stop depth exists and we are in the layer where we have to stop
        if (stop_depth_m and
            math::between(layer.top_depth, stop_depth_m.value(), layer.bot_depth)) {
            stop_depth_was_reached = true;
            return stop_depth_m.value();
        } else {
            return seismo::ray_direction_down(pz0) ? layer.bot_depth : layer.top_depth;
        }
    }();
    auto s_end = (z_interface - z0) / (c * pz0);
    size_t num_steps = std::floor(s_end / ds);
    auto s_step = s_end / num_steps;
    // s has to start at 0 because this calculation is done starting from the current point of
    // the ray.
    // one element more since first step (for s = 0) should also be stored.
    std::vector<state_type> states_vector;
    std::vector<double> arclengths;
    states_vector.reserve(num_steps + 1);
    arclengths.reserve(num_steps + 1);
    const auto x_end = x0 + s_end * c * px0;
    const auto y_end = y0 + s_end * c * py0;
    if (not model.in_horizontal_extent(x_end, y_end)) {
        // ray left model to the side and did not reach top/bottom of layer

        return {Status::OutOfBounds, {}};
    }
    for (auto i = 0U; i < num_steps + 1; ++i) {
        auto arclength_in_layer = i * s_step;
        arclengths.emplace_back(s_start + arclength_in_layer);
        auto x = x0 + arclength_in_layer * c * px0;
        auto y = y0 + arclength_in_layer * c * py0;
        auto z = z0 + arclength_in_layer * c * pz0;
        states_vector.emplace_back(state_type{x, y, z, px0, py0, pz0, t0 + arclength_in_layer / c});
    }
    // "fix" numerical inaccuracy
    states_vector.back()[Index::Z] = z_interface;
    return {stop_depth_was_reached ? Status::StopDepthReached : Status::Success,
            RaySegment{states_vector, arclengths}};
}


RayTracer::RayTracer(const VelocityModel& velocity_model) : model(std::move(velocity_model)) {}

RayTracingResult<Ray> RayTracer::trace_ray(state_type initial_state,
                                           const std::vector<WaveType>& ray_code,
                                           std::optional<double> stop_depth, double step_size,
                                           double max_step) {
    if (not model.in_model(initial_state[Index::X], initial_state[Index::Y],
                           initial_state[Index::Z])) {
        throw std::domain_error(impl::Formatter()
                                << "Point "
                                << point_to_str(initial_state[Index::X], initial_state[Index::Y],
                                                initial_state[Index::Z])
                                << " not in model " << model);
    }
    stop_depth_m = stop_depth;
    auto layer_index = model.layer_index(initial_state[Index::Z]).value();
    auto current_layer = model[layer_index];
    auto segment = trace_layer(initial_state, current_layer, 0., step_size, max_step);
    if (segment.status == Status::OutOfBounds) {
        // early exit when ray left model horizontally
        return {Status::OutOfBounds, {}};
    }
    if (segment.status == Status::StopDepthReached) {
        return {Status::StopDepthReached, Ray{{segment.value()}}};
    }
    if (ray_code.empty()) {
        // only trace one layer for empty ray code
        return {Status::Success, {{{segment.value()}}}};
    }
    Ray ray{{segment.value()}};
    for (auto ray_type : ray_code) {
        // TODO this line and the trace_layer call below violates the law of Demeter, try to
        //  refactor it by improving Ray class
        auto [x, y, z, px, py, pz, t] = ray.segments.back().data.back();
        // reflected waves stay in the same layer, so the index doesn't change
        if (ray_type == WaveType::Transmitted) {
            layer_index += seismo::ray_direction_down(pz) ? 1 : -1;
            current_layer = model[layer_index];
        }
        if (layer_index < 0 or layer_index >= model.size()) {
            throw std::runtime_error(impl::Formatter() << "Ray left model at top or bottom.");
        }
        auto [v_above, v_below] = model.interface_velocities(z);
        auto [px_new, py_new, pz_new] = snells_law(px, py, pz, v_above, v_below, ray_type);
        state_type new_initial_state{x, y, z, px_new, py_new, pz_new, t};
        segment = trace_layer(new_initial_state, current_layer,
                              ray.segments.back().arclength.back(), step_size, max_step);
        if (segment.status == Status::OutOfBounds) {
            return {Status::OutOfBounds, {}};
        }
        ray.segments.push_back(segment.value());
        if (segment.status == Status::StopDepthReached) {
            return {Status::StopDepthReached, ray};
        }
    }
    return {Status::Success, ray};
}


RayTracingResult<RaySegment> RayTracer::trace_layer(const state_type& initial_state,
                                                    const Layer& layer, double s_start, double ds,
                                                    double max_ds) {
    if (layer.gradient == 0) {
        return trace_layer_const(initial_state, layer, s_start, ds);
    } else {
        return trace_layer_gradient(initial_state, layer, s_start, ds, max_ds);
    }
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
    std::pair<Eigen::Tensor2cd, Eigen::Tensor2cd>
    transform(const Eigen::Matrix2cd& P, const Eigen::Matrix2cd& Q, WaveType wave_type,
              const state_type& old_state, const state_type& new_state, int layer_index,
              const VelocityModel& model) {
        // TODO modify interface unit vector (params x2, y2, z2) for more general velocity model.
        //  Here it is assumed the model consists only of horizontal layers.
        msg(layer_index);
        msg(wave_type);
        auto i_S =
            math::angle(old_state[Index::PX], old_state[Index::PY], old_state[Index::PZ], 0, 0, 1);
        msg(i_S);
        auto i_R = wave_type == WaveType::Transmitted
                       ? math::angle(new_state[Index::PX], new_state[Index::PY],
                                     new_state[Index::PZ], 0, 0, 1)
                       : i_S;
        msg(i_R);
        // epsilon is introduced by eq. 2.4.71, Cerveny2001. This formula is simplified for
        // horizontal interfaces (unit vector (0, 0, 1)).
        auto epsilon = std::copysign(1., old_state[Index::PZ]);
        msg(epsilon);
        // for a downgoing transmitted ray the velocity above the interface is the before
        // velocity and the velocity below the interface is the after velocity.
        auto [V_top, V_bottom] = model.interface_velocities(old_state[Index::Z]);
        auto V_before = V_top, V_after = V_bottom;
        if (wave_type == WaveType::Reflected) {
            V_after = V_before;
        } else {
            if (not seismo::ray_direction_down(old_state[Index::PZ])) {
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
        auto old_gradient = model[layer_index].gradient;
        msg(old_gradient);
        auto next_layer_index =
            seismo::next_layer_index(layer_index, old_state[Index::PZ], wave_type);
        auto new_gradient = model[next_layer_index].gradient;
        msg(new_gradient);
        // eq. (4.4.53) from Cerveny2001
        auto E = E_(V_before, i_S, epsilon, old_gradient);
        msg(E);
        auto E_tilde = E_tilde_(wave_type, V_after, i_R, epsilon, new_gradient);
        msg(E_tilde);
        auto u = u_(wave_type, V_before, V_after, i_S, i_R, epsilon);
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
        return {Eigen::TensorMap<Eigen::Tensor2cd>(P_tilde.data(), {2, 2}),
                Eigen::TensorMap<Eigen::Tensor2cd>(Q_tilde.data(), {2, 2})};
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

RayTracingResult<Beam> RayTracer::trace_beam(state_type initial_state, Meter beam_width,
                                             AngularFrequency beam_frequency,
                                             const std::vector<WaveType>& ray_code,
                                             std::optional<double> stop_depth, double step_size,
                                             double max_step) {
    // first trace ray kinematically
    auto ray = trace_ray(initial_state, ray_code, stop_depth, step_size, max_step);
    if (not ray.result) {
        // ray tracing failed
        return {ray.status, {}};
    }
    InterfacePropagator ip;
    Beam beam(beam_width, beam_frequency);
    auto [x, y, z, px, py, pz, t] = initial_state;
    auto v0 = model.eval_at(x, y, z).value();
    using cdouble = std::complex<double>;
    // initial values for P, Q
    Eigen::Tensor3cd P0(1, 2, 2);
    P0.setValues({{{1j / v0, 0}, {0, 1j / v0}}});
    Eigen::Tensor3cd Q0(1, 2, 2);
    Q0.setValues({{{beam_frequency.get() * beam_width.get() * beam_width.get() / v0, 0},
                   {0, beam_frequency.get() * beam_width.get() * beam_width.get() / v0}}});
    std::ptrdiff_t segment_index = 0;
    for (const auto& segment : ray.value()) {
        auto [x, y, z, px, py, pz, t] = segment.data.front();
        auto layer_index = model.layer_index(x, y, z).value();
        // evaluate velocity at all points of the ray
        std::vector<double> v;
        v.reserve(segment.data.size());
        std::transform(
            segment.data.begin(), segment.data.end(), std::back_inserter(v),
            [&](const state_type& state) {
                return model.eval_at(state[Index::X], state[Index::Y], state[Index::Z]).value();
            });
        auto sigma_ = math::cumtrapz(v, segment.arclength, 0.);
        Eigen::TensorRef<Eigen::Tensor1d> sigma =
            Eigen::TensorMap<Eigen::Tensor1d>(sigma_.data(), sigma_.size());
        // TODO this could be optimized since P0 is constant in my use case to store it only once
        Eigen::Tensor3cd P(sigma_.size(), 2, 2);
        P = P0.broadcast(std::array<size_t, 3>{sigma_.size(), 1, 1});
        Eigen::Tensor3cd Q(sigma_.size(), 2, 2);
        Q = Q0.broadcast(std::array<size_t, 3>{sigma_.size(), 1, 1}) +
            sigma.cast<cdouble>()
                    .reshape(std::array<size_t, 3>{sigma_.size(), 1, 1})
                    .broadcast(std::array<int, 3>{1, 2, 2}) *
                P0.broadcast(std::array<size_t, 3>{sigma_.size(), 1, 1});
        beam.segments.emplace_back(segment, P, Q, v);
        if (segment_index < (ray.value().size() - 1)) {
            // if we are not at the last segment of the ray, transform dynamic ray tracing across
            // interface (calculate new P0, Q0)
            auto wave_type = ray_code[segment_index];
            auto new_initial_state = ray.value()[segment_index + 1].data.front();
            auto [P0_new, Q0_new] =
                ip.transform(last_element(P), last_element(Q), wave_type, segment.data.back(),
                             new_initial_state, layer_index, model);
            P0.reshape(std::array<long, 2>{2, 2}) = P0_new;
            Q0.reshape(std::array<long, 2>{2, 2}) = Q0_new;
        }
        ++segment_index;
    }
    return {Status::Success, beam};
}

RayTracingResult<Ray> RayTracer::trace_ray(state_type initial_state, const std::string& ray_code,
                                           std::optional<double> stop_depth, double step_size,
                                           double max_step) {
    return trace_ray(initial_state, seismo::make_ray_code(ray_code), stop_depth, step_size,
                     max_step);
}

RayTracingResult<Beam> RayTracer::trace_beam(state_type initial_state, Meter beam_width,
                                             AngularFrequency beam_frequency,
                                             const std::string& ray_code,
                                             std::optional<double> stop_depth, double step_size,
                                             double max_step) {
    return trace_beam(initial_state, beam_width, beam_frequency, seismo::make_ray_code(ray_code),
                      stop_depth, step_size, max_step);
}
