#include "twopoint.hpp"
#include "utils.hpp"
#include <cmath>
#include <iostream>
#include <valarray>

#define USEDEBUG false
#if USEDEBUG
#define msg(x) std::cout << #x << " = " << x << std::endl;
#else
#define msg(x)
#endif

/**
 * Summ values in array treating nan as zero.
 * @param in
 * @return
 */
TwoPointRayTracing::array_t::value_type nansum(const TwoPointRayTracing::array_t& in) {
    TwoPointRayTracing::array_t::value_type sum{0};
    for (auto e : in) {
        if (std::isnan(e)) {
            break;
        }
        sum += e;
    }
    return sum;
}


TwoPointRayTracing::TwoPointRayTracing(const VelocityModel& velocity_model) :
        model(velocity_model),
        num_layers(model.size()),
        a(model.size() + 1),
        b(model.size() + 1),
        z(model.interface_depths().data(), model.interface_depths().size()) {
    // index zero is reserved for property of source layer, will be filled in trace method
    size_t index = 1;
    for (const auto& layer : model) {
        a[index] = layer.gradient;
        b[index] = layer.intercept;
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
template <typename A>
TwoPointRayTracing::array_t X_tilde(double q, const A& delta_a, const A& mu_tilde,
                                    const A& epsilon_tilde, const A& omega_tilde,
                                    const A& h_tilde) {
    return delta_a * mu_tilde *
               (std::sqrt(std::pow(q, -2) + epsilon_tilde) -
                std::sqrt(std::pow(q, -2) + omega_tilde)) +
           (1. - delta_a) * h_tilde / std::sqrt(std::pow(q, -2) + epsilon_tilde);
}

// eq. B9
template <typename A>
TwoPointRayTracing::array_t X_tilde_prime(double q, const A& delta_a, const A& mu_tilde,
                                          const A& epsilon_tilde, const A& omega_tilde,
                                          const A& h_tilde) {
    return delta_a * mu_tilde / (q * q) *
               (1 / std::sqrt(1 + omega_tilde * q * q) - 1 / std::sqrt(1 + epsilon_tilde * q * q)) +
           (1. - delta_a) * h_tilde / std::pow(1 + epsilon_tilde * q * q, 1.5);
}

// eq. B10
template <typename A>
TwoPointRayTracing::array_t X_tilde_double_prime(double q, const A& delta_a, const A& mu_tilde,
                                                 const A& epsilon_tilde, const A& omega_tilde,
                                                 const A& h_tilde) {
    return delta_a * mu_tilde / (q * q * q) *
               ((2 + 3. * epsilon_tilde * q * q) / std::pow(1 + epsilon_tilde * q * q, 1.5) -
                (2 + 3. * omega_tilde * q * q) / std::pow(1 + omega_tilde * q * q, 1.5)) -
           (1. - delta_a) * 3 * h_tilde * epsilon_tilde * q /
               std::pow(1 + epsilon_tilde * q * q, 2.5);
}

// eq. 16
template <typename A>
double f_tilde(double q, double X, const A& delta_a, const A& mu_tilde, const A& epsilon_tilde,
               const A& omega_tilde, const A& h_tilde) {
    return nansum(X_tilde(q, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde)) - X;
}

// eq. B3
template <typename A>
double f_tilde_prime(double q, const A& delta_a, const A& mu_tilde, const A& epsilon_tilde,
                     const A& omega_tilde, const A& h_tilde) {
    return nansum(X_tilde_prime(q, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde));
}

// eq. B4
template <typename A>
double f_tilde_double_prime(double q, const A& delta_a, const A& mu_tilde, const A& epsilon_tilde,
                            const A& omega_tilde, const A& h_tilde) {
    return nansum(X_tilde_double_prime(q, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde));
}

template <typename AA>
double next_q(double q, double X, const AA& delta_a, const AA& mu_tilde, const AA& epsilon_tilde,
              const AA& omega_tilde, const AA& h_tilde) {
    double A =
        0.5 * f_tilde_double_prime(q, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde);
    double B = f_tilde_prime(q, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde);
    double C = f_tilde(q, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde);
    double delta_q_plus = (-B + std::sqrt(B * B - 4 * A * C)) / (2 * A);
    double delta_q_minus = (-B - std::sqrt(B * B - 4 * A * C)) / (2 * A);
    // both q plus and q minus are 0D tensors, get their value out to pass to function
    double q_plus = (q + delta_q_plus);
    double q_minus = (q + delta_q_minus);
    // use all to convert 0D array of bool to bool explicitly
    if (std::abs(f_tilde(q_plus, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde)) <
        std::abs(f_tilde(q_minus, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde))) {
        return q_plus;
    } else {
        return q_minus;
    }
}

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
TwoPointRayTracing::array_t guard(const TwoPointRayTracing::array_t& e,
                                  const TwoPointRayTracing::array_t& guard) {
    TwoPointRayTracing::array_t out(e.size());
    for (size_t i = 0; i < e.size(); ++i) {
        out[i] = guard[i] == 0 ? guard[i] : e[i];
    }
    return out;
}


double q_to_p(double q, double vM) {
    return std::sqrt(q * q / (vM * vM + vM * vM * q * q));
}


TwoPointRayTracing::array_t delta(const TwoPointRayTracing::array_t& in) {
    return in.apply([](double d) -> double { return d != 0 ? 1 : 0; });
}

slowness_t TwoPointRayTracing::trace(position_t source, position_t receiver, double accuracy) {
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
    auto source_index = static_cast<std::size_t>(model.layer_index(source_z).value() + 1);
    if (source_index > num_layers) {
        source_index = num_layers;
    }
    // insert source layer properties in first place
    a[0] = a[source_index];
    b[0] = b[source_index];

    // eq. A5
    // TODO move instantiation to constructor so they can be reused.
    auto epsilon = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        epsilon[k] = std::pow(a[k] * z[k - 1] + b[k], 2);
    }
    epsilon[0] = std::pow(a[source_index] * z[source_index - 1] + b[source_index], 2);

    // eq. A6
    // TODO decide if kept as loop or as xtensor view initialization as for epsilon
    auto omega = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        omega[k] = std::pow(a[k] * z[k] + b[k], 2);
    }
    omega[0] = std::pow(a[source_index] * source_z + b[source_index], 2);

    // eq. A7
    auto h = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        h[k] = (a[k] * z[k - 1] + b[k]) * (z[k] - z[k - 1]);
    }
    h[0] = (a[source_index] * z[source_index - 1] + b[source_index]) *
           (source_z - z[source_index - 1]);

    bool source_below_receiver = source_z > receiver_z;
    auto mu = array_t(num_layers + 1);
    for (size_t k = 0; k <= num_layers; k++) {
        mu[k] = mu_k(k, source_index, num_layers, source_below_receiver);
    }

    msg("Initialized");

    double vM = highest_velocity_between(source_z, receiver_z, model);

    // eq. A10
    array_t mu_tilde = mu * vM / a;
    for (size_t i = 0; i < mu_tilde.size(); ++i) {
        if (a[i] == 0) {
            mu_tilde[i] = 0;
        }
    }

    // eq. A11
    array_t h_tilde = mu * h / vM;


    // eq. A12
    array_t epsilon_tilde = 1 - epsilon / (vM * vM);


    // eq. A13
    array_t omega_tilde = 1 - omega / (vM * vM);


    // eq. A9
    array_t delta_a = delta(a);


    // eq. C13
    auto d1 =
        nansum(delta_a * 0.5 * mu_tilde * (epsilon_tilde - omega_tilde) + (1. - delta_a) * h_tilde);

    // eq. C18
    array_t delta_epsilon = delta(epsilon_tilde);

    // eq. C19
    array_t delta_omega = delta(omega_tilde);

    // eq. C14
    auto c0 = nansum(
        delta_a * mu_tilde *
            (delta_epsilon * std::sqrt(epsilon_tilde) - delta_omega * std::sqrt(omega_tilde)) +
        (1. - delta_a) * delta_epsilon * h_tilde / std::sqrt(epsilon_tilde));


    // eq. C16
    auto cminus1 = nansum(delta_a * (delta_omega - delta_epsilon) * mu_tilde);


    // eq. C17
    auto cminus2 =
        nansum(delta_a * 0.5 * mu_tilde *
                   (guard(delta_epsilon / std::sqrt(epsilon_tilde), delta_epsilon) -
                    guard(delta_omega / std::sqrt(omega_tilde), delta_omega)) -
               (1. - delta_a) * delta_epsilon * h_tilde * 0.5 / std::pow(epsilon_tilde, 1.5));

    // horizontal distance between source and receiver
    auto X = std::sqrt(std::pow(source_x - receiver_x, 2) + std::pow(source_y - receiver_y, 2));
    msg(X);
    auto alpha1 = d1;
    msg(alpha1);
    auto alpha2 = c0 * (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    msg(alpha2);
    auto beta1 = (c0 * cminus1 + d1 * cminus2) / (c0 * cminus2 - cminus1 * cminus1);
    msg(beta1);
    auto beta2 = (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    msg(beta2);
    auto numerator = (beta1 * X - alpha1 +
                      sqrt((beta1 * beta1 - 4 * beta2) * X * X +
                           2 * (2 * alpha2 - alpha1 * beta1) * X + alpha1 * alpha1));
    msg(numerator);
    auto denominator = 2 * (alpha2 - beta2 * X);
    msg(denominator);
    double q = numerator / denominator;
    msg(q);
    double q_next;
    msg("Before loop");
    while (std::isfinite(q) and q != 0) {
        q_next = next_q(q, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde);
        msg(q_next);
        if (std::abs(q - q_next) < accuracy) {
            q = q_next;
            break;
        }
        q = q_next;
    }
    double horizontal_slowness = q_to_p(q, vM);
    double c = model.eval_at(source_x, source_y, source_z).value();
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