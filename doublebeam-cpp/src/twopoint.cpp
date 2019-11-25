#include "twopoint.hpp"
#include "utils.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <valarray>


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& array) {
    os << "(";
    for (auto item : array) {
        os << item << " ";
    }
    os << ")";
    return os;
}

#define USEDEBUG true
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
            continue;
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
    TwoPointRayTracing::array_t tmp1 = mu_tilde / q *
                                       (std::sqrt(1 + std::pow(q, 2) * epsilon_tilde) -
                                        std::sqrt(1 + std::pow(q, 2) * omega_tilde));
    TwoPointRayTracing::array_t tmp2 = h_tilde * q / std::sqrt(1 + std::pow(q, 2) * epsilon_tilde);
    // this simulates multiplication by delta_a but works for Nan
    tmp1[delta_a == 0] = 0;
    // simulate multiplication by 1 - delta_a
    tmp2[delta_a == 1] = 0;
    return tmp1 + tmp2;
}

// eq. B9
template <typename A>
TwoPointRayTracing::array_t X_tilde_prime(double q, const A& delta_a, const A& mu_tilde,
                                          const A& epsilon_tilde, const A& omega_tilde,
                                          const A& h_tilde) {
    TwoPointRayTracing::array_t tmp1 =
        mu_tilde / (q * q) *
        (1 / std::sqrt(1 + omega_tilde * q * q) - 1 / std::sqrt(1 + epsilon_tilde * q * q));
    TwoPointRayTracing::array_t tmp2 = h_tilde / std::pow(1 + epsilon_tilde * q * q, 1.5);
    tmp1[delta_a == 0] = 0;
    tmp2[delta_a == 1] = 0;
    return tmp1 + tmp2;
}

// eq. B10
template <typename A>
TwoPointRayTracing::array_t X_tilde_double_prime(double q, const A& delta_a, const A& mu_tilde,
                                                 const A& epsilon_tilde, const A& omega_tilde,
                                                 const A& h_tilde) {
    TwoPointRayTracing::array_t tmp1 =
        mu_tilde / (q * q * q) *
        ((2 + 3. * epsilon_tilde * q * q) / std::pow(1 + epsilon_tilde * q * q, 1.5) -
         (2 + 3. * omega_tilde * q * q) / std::pow(1 + omega_tilde * q * q, 1.5));
    TwoPointRayTracing::array_t tmp2 =
        3 * h_tilde * epsilon_tilde * q / std::pow(1 + epsilon_tilde * q * q, 2.5);
    tmp1[delta_a == 0] = 0;
    tmp2[delta_a == 1] = 0;
    return tmp1 - tmp2;
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

// Solve quadratic equation (7) using eq. (11) only for q instead of p.
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
    if (std::abs(f_tilde(q_plus, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde)) <
        std::abs(f_tilde(q_minus, X, delta_a, mu_tilde, epsilon_tilde, omega_tilde, h_tilde))) {
        return q_plus;
    } else {
        return q_minus;
    }
}

// Transformed version of eq. (15) for p
double q_to_p(double q, double vM) {
    return std::sqrt(q * q / (vM * vM + vM * vM * q * q));
}


template <typename T>
T Heaviside(T x) {
    return x < 0 ? T{0.} : T{1};
}


// Calculate starting value for q using either C1 or C6 depending on value of c1
double initial_q(double X, double d1, double c0, double c1, double cminus1, double cminus2) {
    auto alpha1 = d1;
    msg(alpha1);
    auto alpha2 = c0 * (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    msg(alpha2);
    auto beta1 = (c0 * cminus1 + d1 * cminus2) / (c0 * cminus2 - cminus1 * cminus1);
    msg(beta1);
    auto beta2 = (c0 * c0 + d1 * cminus1) / (cminus1 * cminus1 - c0 * cminus2);
    msg(beta2);

    if (c1 == 0) {
        // Use equation C1 to estimate initial value for q
        auto numerator = (beta1 * X - alpha1 +
                          sqrt((beta1 * beta1 - 4 * beta2) * X * X +
                               2 * (2 * alpha2 - alpha1 * beta1) * X + alpha1 * alpha1));
        msg(numerator);
        auto denominator = 2 * (alpha2 - beta2 * X);
        msg(denominator);
        double q = numerator / denominator;
        return q;
    } else {
        // This uses equation C6 to estimate initial value for q. The paper mistakenly states this
        // equation should be used when c_1 = 0, but it should be used when c1 != 0. This was
        // confirmed by email from the corresponding author.
        auto q = X / alpha1 +
                 (X / c1 - X / alpha1 - c0 / c1) * Heaviside(X - (alpha1 * c0) / (alpha1 - c1));
        return q;
    }
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
    msg(source_index);
    // insert source layer properties in first place
    a[0] = a[source_index];
    b[0] = b[source_index];
    msg(a);
    msg(b);
    msg(z);
    // eq. A5
    // TODO Many arrays only depend on the source depth. For my use case, this will be constant so
    //  they could be precalculated.
    auto epsilon = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        epsilon[k] = std::pow(a[k] * z[k - 1] + b[k], 2);
    }
    epsilon[0] = std::pow(a[source_index] * z[source_index - 1] + b[0], 2);
    msg(epsilon);
    // eq. A6
    auto omega = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        omega[k] = std::pow(a[k] * z[k] + b[k], 2);
    }
    omega[0] = std::pow(a[source_index] * source_z + b[source_index], 2);
    msg(omega);
    // eq. A7
    auto h = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        h[k] = (a[k] * z[k - 1] + b[k]) * (z[k] - z[k - 1]);
    }
    h[0] = (a[source_index] * z[source_index - 1] + b[source_index]) *
           (source_z - z[source_index - 1]);
    msg(h);
    bool source_below_receiver = source_z > receiver_z;
    auto mu = array_t(num_layers + 1);
    for (size_t k = 0; k <= num_layers; k++) {
        mu[k] = mu_k(k, source_index, num_layers, source_below_receiver);
    }
    msg(mu);

    double vM = highest_velocity_between(source_z, receiver_z, model);
    msg(vM);
    // eq. A10
    array_t mu_tilde = mu * vM / a;
    // mu_tilde is only ever used multiplied together with delta_a as a guard, so set the relevant
    // values to zero early to deal with inf.
    // This also means mu_tilde has delta_a already multiplied in hereafter
    mu_tilde[a == 0] = 0;
    msg(mu_tilde);
    // eq. A11
    array_t h_tilde = mu * h / vM;
    msg(h_tilde);

    // eq. A12
    array_t epsilon_tilde = 1 - epsilon / (vM * vM);
    msg(epsilon_tilde);

    // eq. A13
    array_t omega_tilde = 1 - omega / (vM * vM);
    msg(omega_tilde);

    // eq. A9
    array_t delta_a = delta(a);
    msg(delta_a);

    // eq. C13
    array_t d1_unsummed = 0.5 * mu_tilde * (epsilon_tilde - omega_tilde) + (1. - delta_a) * h_tilde;
    msg(d1_unsummed);
    auto d1 = nansum(d1_unsummed);
    msg(d1);
    // eq. C18
    array_t delta_epsilon = delta(epsilon_tilde);
    msg(delta_epsilon);
    // eq. C19
    array_t delta_omega = delta(omega_tilde);
    msg(delta_omega);
    //    array_t tmp1 = sqrt(epsilon_tilde);
    //    tmp1[delta_epsilon == 0] = 0;
    //    array_t tmp2 = sqrt(omega_tilde);
    //    tmp2[delta_omega == 0] = 0;
    //    array_t tmp3 = mu_tilde * (tmp1 - tmp2);
    //    tmp3[delta_a == 0] = 0;
    // eq. C14
    // TODO fix delta/guards here
    array_t c0_unsummed =
        delta_a * mu_tilde *
            (delta_epsilon * std::sqrt(epsilon_tilde) - delta_omega * std::sqrt(omega_tilde)) +
        (1. - delta_a) * delta_epsilon * h_tilde / std::sqrt(epsilon_tilde);
    msg(c0_unsummed);
    auto c0 = nansum(c0_unsummed);
    msg(c0);

    // eq. C16
    array_t cminus1_unsummed = delta_a * (delta_omega - delta_epsilon) * mu_tilde;
    msg(cminus1_unsummed);
    auto cminus1 = nansum(cminus1_unsummed);
    msg(cminus1);

    // eq. C17
    array_t cminus2_unsummed =
        delta_a * 0.5 * mu_tilde *
            (delta_epsilon / std::sqrt(epsilon_tilde) - delta_omega / std::sqrt(omega_tilde)) -
        (1. - delta_a) * delta_epsilon * h_tilde * 0.5 / std::pow(epsilon_tilde, 1.5);
    msg(cminus2_unsummed);
    auto cminus2 = nansum(cminus2_unsummed);
    msg(cminus2);
    // horizontal distance between source and receiver
    auto X = std::sqrt(std::pow(source_x - receiver_x, 2) + std::pow(source_y - receiver_y, 2));
    msg(X);

    // eq. C15
    auto c1 = nansum((1 - delta_a) * (1 - delta_epsilon) * h_tilde);
    msg(c1);

    auto q = initial_q(X, d1, c0, c1, cminus1, cminus2);
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