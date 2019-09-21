#include <complex>

#include "xtensor/xbroadcast.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>

#include "dynamic_raytracing.hpp"
#include "utils.hpp"


DynamicRayTracer::DynamicRayTracer(const VelocityModel& model) : kinematic(model), model(model) {}

using complex = std::complex<double>;


class InterfacePropagator {

    using matrix_t = xt::xtensor<std::complex<double>, 3>;

public:
    std::pair<matrix_t, matrix_t> transform(matrix_t P, matrix_t Q, char wave_type,
                                            const state_type& old_state,
                                            const state_type& new_state, int layer_index,
                                            const VelocityModel& model) {
        // TODO modify interface unit vector (params x2, y2, z2) for more general velocity model.
        //  Here it is assumed the model consists only of horizontal layers.
        auto i_S =
            math::angle(old_state[Index::PX], old_state[Index::PY], old_state[Index::PZ], 0, 0, 1);
        auto i_R =
            math::angle(new_state[Index::PX], new_state[Index::PY], new_state[Index::PZ], 0, 0, 1);
        // epsilon is introduced by eq. 2.4.71, Cerveny2001. This formula is simplified for
        // horizontal interfaces (unit vector (0, 0, 1)).
        auto epsilon = std::copysign(1., old_state[Index::PZ]);
        // for a downgoing transmitted ray the velocity above the interface is the before velocity
        // and the velocity below the interface is the after velocity.
        auto [V_top, V_bottom] = model.interface_velocities(old_state[Index::Z]);
        auto V_before = V_top, V_after = V_bottom;
        if (wave_type == 'R') {
            V_after = V_before;
        } else {
            if (not seismo::ray_direction_down(old_state[Index::PZ])) {
                std::swap(V_before, V_after);
            }
        }
        // TODO this kappa is only valid for simple velocity model v = v(z) and horizontal
        //  interfaces
        auto kappa = 0.;
        auto cos_kappa = std::cos(kappa), sin_kappa = std::sin(kappa);
        // right equations of (4.4.49) in Cerveny2001
        matrix_t G_orthogonal{{{cos_kappa, -sin_kappa}, {sin_kappa, cos_kappa}}};
        auto G_orthogonal_tilde = G_orthogonal;
        // left equations of (4.4.49) in Cerveny2001
        matrix_t G_parallel = xt::eye(2);
        G_parallel(0) = epsilon * std::cos(i_S);
        auto G_parallel_tilde = G_parallel;
        G_parallel_tilde(0) *= wave_type == 'T' ? 1 : -1;
        // equation (4.4.48) from Cerveny2001
        auto G = xt::linalg::dot(G_parallel, G_orthogonal);
        auto G_tilde = xt::linalg::dot(G_parallel_tilde, G_orthogonal_tilde);

        auto G_inverted = xt::linalg::inv(G);
        auto G_tilde_inverted = xt::linalg::inv(G_tilde);
        // eq. (4.4.53) from Cerveny2001
        double old_gradient = 0, new_gradient = 0;
        auto E = E_(V_before, i_S, epsilon, old_gradient);
        auto E_tilde = E_tilde_(wave_type, V_after, i_R, epsilon, new_gradient);
        auto u = u_(wave_type, V_before, V_after, i_R, epsilon, new_gradient);
        auto D = D_();
        // eq. (4.4.67) Cerveny2001
        namespace xtl = xt::linalg;
        auto P_tilde = xtl::dot(
            G_tilde_inverted,
            xtl::dot(G, P) + xtl::dot(E - E_tilde - u * D, xtl::dot(xt::transpose(G_inverted), Q)));
        auto Q_tilde = xtl::dot(xt::transpose(G_tilde), xtl::dot(xt::transpose(G_inverted), Q));
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
        return {{{E11, E12}, {E12, E22}}};
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
    matrix_t E_tilde_(char wave_type, double V_tilde, double i_R, double epsilon,
                      double new_gradient) const {
        auto dV_tilde_dz1 = 0.;
        auto dV_tilde_dz2 = 0.;
        auto dV_tilde_dz3 = new_gradient;
        auto minus_plus = wave_type == 'R' ? -1. : 1.;
        auto E11 = -sin(i_R) / (V_tilde * V_tilde);
        E11 *= ((1 + std::pow(cos(i_R), 2)) * dV_tilde_dz1 +
                minus_plus * epsilon * std::cos(i_R) * std::sin(i_R) * dV_tilde_dz3);
        auto E12 = -sin(i_R) / (V_tilde * V_tilde) * dV_tilde_dz2;
        auto E22 = 0.;
        return {{{E11, E12}, {E12, E22}}};
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
    static double u_(char wave_type, double V, double V_tilde, double i_S, double i_R,
                     double epsilon) {
        auto minusplus = wave_type == 'T' ? -1 : 1;
        return epsilon * (std::cos(i_S) / V + minusplus * std::cos(i_R) / V_tilde);
    }

    /**
     * Eq. 4.4.15 from Cerveny2001
     * For the currently implemented velocity layer with horizontal interfaces only, this function
     * is zero everywhere since it contains the second derivative of the interface function Sigma
     * in the numerator. Sigma = Sigma(z3) for horizontal interfaces.
     */
    static matrix_t D_() {
        return xt::zeros<complex>({2, 2});
    }
};


Beam DynamicRayTracer::trace_beam(state_type initial_state, double beam_width,
                                  double beam_frequency, const std::string& ray_code,
                                  double step_size, double max_step) {
    if (not model.in_model(initial_state[Index::Z])) {
        std::stringstream s;
        auto [top, bottom] = model.get_top_bottom();
        s << "Starting point of beam (" << initial_state[Index::X] << " " << initial_state[Index::Y]
          << " " << initial_state[Index::Z] << ") outside of model ( " << top << " " << bottom
          << ").";
        throw std::domain_error(s.str());
    }
    Beam beam(beam_width, beam_frequency);
    InterfacePropagator ip;
    // initial values for P, Q
    auto v0 = model.eval_at(initial_state[Index::Z]);
    xt::xtensor<complex, 3> P0{{{1j / v0, 0}, {0, 1j / v0}}};
    xt::xtensor<complex, 3> Q0{{{beam_frequency * beam_width * beam_width / v0, 0},
                                {0, beam_frequency * beam_width * beam_width / v0}}};
    // trace ray through model to get the path and slowness
    auto ray = kinematic.trace_ray(initial_state, ray_code, step_size, max_step);
    auto layer_index = model.layer_index(initial_state[Index::Z]);
    auto layer_indices =
        seismo::ray_code_to_layer_indices(ray_code, initial_state[Index::PZ], layer_index);
    for (size_t segment_index = 0; segment_index < ray.size(); ++segment_index) {
        auto segment = ray[segment_index];
        layer_index = layer_indices[segment_index];
        // evaluate velocity at all points of the ray
        // TODO this was already done during ray tracing itself, maybe we can reuse the results
        std::vector<double> v;
        //v.reserve(segment.data.size());
        std::transform(segment.data.begin(), segment.data.end(), std::back_inserter(v),
                       [&](const state_type& state) { return model.eval_at(state[Index::Z]); });
        auto sigma_ = math::cumtrapz(v.begin(), v.end(), segment.arclength.begin(),
                                     segment.arclength.end(), 0.);
        auto sigma = xt::adapt(sigma_, {sigma_.size(), 1UL, 1UL});
        xt::xtensor<complex, 3> Q = Q0 + sigma * P0;
        xt::xtensor<complex, 3> P = xt::broadcast(P0, {sigma.size(), 2UL, 2UL});
        std::cout << "Q: " << Q << std::endl;
        std::cout << "P: " << P << std::endl;
        beam.segments.emplace_back(segment, P, Q);
        // set P0, Q0 to new initial state after crossing the interface
        // std::bind(P0, Q0) = ip.transform(P0, xt::view(Q, xt::keep(-1)), );
    }
    return beam;
}
