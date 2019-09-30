#include "twopoint.hpp"
#include "utils.hpp"
#include <cmath>
#include <iostream>
#include <xtensor/xadapt.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xio.hpp>


TwoPointRayTracing::TwoPointRayTracing(VelocityModel& velocity_model) :
        model(velocity_model),
        num_layers(model.size()),
        gradients(model.size() + 1),
        intercepts(model.size() + 1),
        interface_depths(model.interface_depths()) {
    // index zero is reserved for property of source layer, will be filled in trace method
    size_t index = 1;
    for (const auto& layer : model) {
        gradients[index] = layer.gradient;
        intercepts[index] = layer.intercept;
        ++index;
    }
}

std::string stringify(position_t pos) {
    auto [x, y, z] = pos;
    return "(" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ")";
}

/**
 * Eq. A8
 * @param k Index of current layer.
 * @param s Index of source layer.
 * @param n Number of layers.
 * @param source_below_receiver Fang2019s algorithm is only designed for
 * sources below the receiver. If the source is above the receiver, this
 * function has to be modified. Mu_k normally returns 1 for all layers above
 * the source layer and 0 for all layers below. This must be reversed if the
 * source is above the receiver.
 * @return
 */
double mu_k(size_t k, size_t s, size_t n, bool source_below_receiver) {
    if (1 <= k and k <= s - 1) {
        return source_below_receiver ? 1 : 0;
    }
    if (s <= k and k <= n) {
        return source_below_receiver ? 0 : 1;
    }
    if (k == 0) {
        return 1 - mu_k(s, s, n, source_below_receiver);
    }
    return 0;
}

// eq. A3
auto TwoPointRayTracing::X_tilde(double q) {
    return delta_a * mu_tilde *
               (xt::sqrt(std::pow(q, -2) + epsilon_tilde) -
                xt::sqrt(std::pow(q, -2) + omega_tilde)) +
           (1 - delta_a) * h_tilde / xt::sqrt(std::pow(q, -2) + epsilon_tilde);
}

// eq. B9
auto TwoPointRayTracing::X_tilde_prime(double q) {
    return delta_a * mu_tilde / (q * q) *
               (1 / xt::sqrt(1 + omega_tilde * q * q) - 1 / xt::sqrt(1 + epsilon_tilde * q * q)) +
           (1 - delta_a) * h_tilde / xt::pow(1 + epsilon_tilde * q * q, 1.5);
}

// eq. B10
auto TwoPointRayTracing::X_tilde_double_prime(double q) {
    return delta_a * mu_tilde / (q * q * q) *
               ((2 + 3 * epsilon_tilde * q * q) / xt::pow(1 + epsilon_tilde * q * q, 1.5) -
                (2 + 3 * omega_tilde * q * q) / xt::pow(1 + omega_tilde * q * q, 1.5)) -
           (1 - delta_a) * 3 * h_tilde * epsilon_tilde * q /
               xt::pow(1 + epsilon_tilde * q * q, 2.5);
}

// eq. 16
double TwoPointRayTracing::f_tilde(double q, double X) {
    return xt::sum(X_tilde(q))[0] - X;
}

// eq. B3
double TwoPointRayTracing::f_tilde_prime(double q) {
    return xt::sum(X_tilde_prime(q))[0];
}

// eq. B4
double TwoPointRayTracing::f_tilde_double_prime(double q) {
    return xt::sum(X_tilde_double_prime(q))[0];
}

double TwoPointRayTracing::next_q(double q, double X) {
    auto A = 0.5 * f_tilde_double_prime(q);
    auto B = f_tilde_prime(q);
    auto C = f_tilde(q, X);
    double delta_q_plus = (-B + std::sqrt(B * B - 4 * A * C)) / (2 * A);
    double delta_q_minus = (-B - std::sqrt(B * B - 4 * A * C)) / (2 * A);
    // both q plus and q minus are 0D tensors, get their value out to pass to function
    double q_plus = (q + delta_q_plus);
    double q_minus = (q + delta_q_minus);
    // use all to convert 0D array of bool to bool explicitly
    if (std::abs(f_tilde(q_plus, X)) < std::abs(f_tilde(q_minus, X))) {
        return q_plus;
    } else {
        return q_minus;
    }
}

template <typename T>
auto is_invalid(T x) {
    return xt::isinf(x) || xt::isnan(x);
}

template <typename E>
void filter_invalid(E& e, typename E::value_type replace = 0.) {
    xt::filtration(e, is_invalid(e)) = replace;
}

template <typename E1, typename E2>
/**
 * Replace parts of expression by 0 as indicated by the zeros in the guard clause.
 * This function is used because the algorithm described in Fang2019 can yield invalid mathematical
 * operations (division by zero for mu_tilde when the velocity of a layer is constant). The authors
 * use "guard clauses", prefixed with delta_*, which are zero at the index of the offending value
 * and 1 elsewhere.
 * If multiplied with the invalid operation, the result of the operation is canceled out.
 * This behaviour is not valid C++, since a Nan or Inf multiplied by zero is not zero but Nan. This
 * problem can be partially solved using nansum, which sums the sequence of values while treating
 * Nans as zero. This relies implicitly on the assumption that every nan would have a guard in the
 * formula. But due to the propagation behaviour of Nans, this is not enough. Every operatoion
 * involving a Nan results in a Nan, such as adding a value to it. In a formula (for one array
 * element)  such as:
 * guard_clause * Nan + valid_value
 * this result in the Nan annulating the valid value, even when the whole array is nansummed.
 * The authors intent this formula to result in valid value.
 * To emulate this behaviour, I use this function as a guard everywhere a divide by zero could occur
 * with the appropriate delta_* guard and a later addition could be annulled by the Nan.
 * @tparam E1 xtensor xexpression
 * @tparam E2 xtensor xexpression
 * @param e expression to evaluate
 * @param guard Guard clause, consisting of 1 or 0.
 * @return Elements of e where the guard clause is unequal 0, else 0.
 */
auto guard(E1&& e, E2&& guard) {
    return xt::where(xt::equal(guard, 0), guard, e);
}


double q_to_p(double q, double vM) {
    return std::sqrt(q * q / (vM * vM + vM * vM * q * q));
}


slowness_t TwoPointRayTracing::trace(position_t source, position_t receiver,
                                     __attribute__((unused)) double accuracy) {
    auto [source_x, source_y, source_z] = source;
    auto [receiver_x, receiver_y, receiver_z] = receiver;
    if (not model.in_model(source_x, source_y, source_z)) {
        throw std::domain_error(impl::Formatter()
                                << "Source at " << stringify(source) << " outside of model.");
    }
    if (not model.in_model(receiver_x, receiver_y, receiver_z)) {
        throw std::domain_error(impl::Formatter()
                                << "Receiver at " << stringify(receiver) << "  outside of model.");
    }
    // while the paper uses source based indexing, C++ doesn't.
    auto source_index = model.layer_index(source_z) + 1;
    if (source_index > num_layers) {
        source_index = num_layers;
    }
    // insert source layer properties in first place
    gradients[0] = gradients[source_index];
    intercepts[0] = intercepts[source_index];

    auto a = xt::adapt(gradients);
    auto b = xt::adapt(intercepts);
    auto z = xt::adapt(interface_depths);

    // eq. A5
    // TODO move instantiation to constructor so they can be reused.
    auto epsilon = xt::empty<double>({num_layers + 1});
    for (size_t k = 1; k <= num_layers; k++) {
        epsilon(k) = std::pow(a(k) * z(k - 1) + b(k), 2);
    }
    epsilon(0) = std::pow(a(source_index) * z(source_index - 1) + b(source_index), 2);

    // eq. A6
    // TODO decide if kept as loop or as xtensor view initialization as for epsilon
    auto omega = xt::empty<double>({num_layers + 1});
    for (size_t k = 1; k <= num_layers; k++) {
        omega(k) = std::pow(a(k) * z(k) + b(k), 2);
    }
    omega(0) = std::pow(a(source_index) * source_z + b(source_index), 2);

    // eq. A7
    auto h = xt::empty<double>({num_layers + 1});
    for (size_t k = 1; k <= num_layers; k++) {
        h(k) = (a(k) * z(k - 1) + b(k)) * (z(k) - z(k - 1));
    }
    h(0) = (a(source_index) * z(source_index - 1) + b(source_index)) *
           (source_z - z(source_index - 1));

    bool source_below_receiver = source_z > receiver_z;
    auto mu = xt::empty<double>({num_layers + 1});
    for (size_t k = 0; k <= num_layers; k++) {
        mu(k) = mu_k(k, source_index, num_layers, source_below_receiver);
    }

    auto vM = highest_velocity_between(source_z, receiver_z, model);

    // eq. A10
    mu_tilde = mu * vM / a;
    xt::filtration(mu_tilde, xt::equal(a, 0)) = 0;


    // eq. A11
    h_tilde = mu * h / vM;


    // eq. A12
    epsilon_tilde = 1 - epsilon / (vM * vM);


    // eq. A13
    omega_tilde = 1 - omega / (vM * vM);


    // eq. A9
    delta_a = xt::abs(xt::sign(a));


    // eq. C13
    auto d1 = xt::nansum(delta_a * 0.5 * mu_tilde * (epsilon_tilde - omega_tilde) +
                         (1 - delta_a) * h_tilde);

    // eq. C18
    delta_epsilon = xt::abs(xt::sign(epsilon_tilde));

    // eq. C19
    delta_omega = xt::abs(xt::sign(omega_tilde));

    // eq. C14
    auto c0 = xt::nansum(
        delta_a * mu_tilde *
            (delta_epsilon * xt::sqrt(epsilon_tilde) - delta_omega * xt::sqrt(omega_tilde)) +
        (1 - delta_a) * delta_epsilon * h_tilde / xt::sqrt(epsilon_tilde));


    // eq. C16
    auto cminus1 = xt::nansum(delta_a * (delta_omega - delta_epsilon) * mu_tilde);


    // eq. C17
    auto cminus2 =
        xt::nansum(delta_a * 0.5 * mu_tilde *
                       (guard(delta_epsilon / xt::sqrt(epsilon_tilde), delta_epsilon) -
                        guard(delta_omega / xt::sqrt(omega_tilde), delta_omega)) -
                   (1 - delta_a) * delta_epsilon * h_tilde * 0.5 / xt::pow(epsilon_tilde, 1.5));


    // horizontal distance between source and receiver
    auto X = std::sqrt(std::pow(source_x - receiver_x, 2) + std::pow(source_y - receiver_y, 2));

    auto alpha1 = d1;
    auto alpha2 = c0 * (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    auto beta1 = (c0 * cminus1 + d1 * cminus2) / (c0 * cminus2 - cminus1 * cminus1);
    auto beta2 = (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    auto numerator = (beta1 * X - alpha1 +
                      sqrt((beta1 * beta1 - 4 * beta2) * X * X +
                           2 * (2 * alpha2 - alpha1 * beta1) * X + alpha1 * alpha1));
    auto denominator = 2 * (alpha2 - beta2 * X);
    double q = (numerator / denominator)[0];

    double q_next;
    while (std::isfinite(q) and q != 0) {
        q_next = next_q(q, X);
        if (std::abs(q - q_next) < accuracy) {
            q = q_next;
            break;
        }
        q = q_next;
    }
    double horizontal_slowness = q_to_p(q, vM);
    double c = model.eval_at(source_x, source_y, source_z);
    double vertical_slowness =
        std::sqrt(std::pow(c, -2) - horizontal_slowness * horizontal_slowness);
    if (source_below_receiver) {
        // ray should travel upward from source in this case
        vertical_slowness *= -1;
    }
    // calculate angle to x axis
    double phi = math::angle_clockwise(receiver_x - source_x, receiver_y - source_y, 1., 0.);
    double px = std::cos(phi) * horizontal_slowness;
    double py = std::sin(phi) * horizontal_slowness;
    return {px, py, vertical_slowness};
}