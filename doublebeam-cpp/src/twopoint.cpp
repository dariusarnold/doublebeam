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
        b(model.size() + 1),
        z(model.interface_depths().data(), model.interface_depths().size()) {
    // index zero is reserved for property of source layer, will be filled in trace method
    size_t index = 1;
    for (const auto& layer : model) {
        if (layer.gradient != 0) {
            throw std::invalid_argument(impl::Formatter()
                                        << "Can't initialize two point ray tracer for model "
                                           "containing layer with linear velocity gradient.");
        }
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
TwoPointRayTracing::array_t X_tilde(double q, const A& epsilon_tilde, const A& h_tilde) {
    TwoPointRayTracing::array_t tmp = h_tilde / std::sqrt(1 / (q * q) + epsilon_tilde);
    return tmp;
}

// eq. B9
template <typename A>
TwoPointRayTracing::array_t X_tilde_prime(double q, const A& epsilon_tilde, const A& h_tilde) {
    TwoPointRayTracing::array_t tmp = h_tilde / std::pow(1 + epsilon_tilde * q * q, 1.5);
    return tmp;
}


// eq. 16
template <typename A>
double f_tilde(double q, double X, const A& epsilon_tilde, const A& h_tilde) {
    return nansum(X_tilde(q, epsilon_tilde, h_tilde)) - X;
}

// eq. B3
template <typename A>
double f_tilde_prime(double q, const A& epsilon_tilde, const A& h_tilde) {
    return nansum(X_tilde_prime(q, epsilon_tilde, h_tilde));
}


// Solve nonlinear equation (5) using Newton algorithm.
template <typename AA>
double next_q(double q, double X, const AA& epsilon_tilde, const AA& h_tilde) {
    double f_prime = f_tilde_prime(q, epsilon_tilde, h_tilde);
    double f = f_tilde(q, X, epsilon_tilde, h_tilde);
    auto q_new = q - f / f_prime;
    return q_new;
}

// Transformed version of eq. (15) for p
double q_to_p(double q, double vM) {
    return std::sqrt(q * q / (vM * vM + vM * vM * q * q));
}


template <typename T>
T Heaviside(T x) {
    return x < 0 ? T{0.} : T{1};
}


double safe_divide(double numerator, double denominator) {
    if (numerator == 0) {
        return numerator;
    }
    return numerator / denominator;
}


// Calculate starting value for q using either C1 or C6 depending on value of c1
double initial_q(double X, double d1, double c0, double c1, double cminus1, double cminus2) {
    auto alpha1 = d1;
    msg(alpha1);
    auto alpha2 = safe_divide(c0 * c0 * c0 + c0 * d1 * cminus1, cminus1 * cminus1 - c0 * cminus2);
    msg(alpha2);
    auto beta1 = safe_divide(c0 * cminus1 + d1 * cminus2, c0 * cminus2 - cminus1 * cminus1);
    msg(beta1);
    auto beta2 = safe_divide(c0 * c0 + d1 * cminus1, cminus1 * cminus1 - c0 * cminus2);
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

double horizontal_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
}

/**
 * This function is not the pure implementation of the algorithm in the paper
 * "A fast and robust two-point ray tracing method in layered media with constant or linearly
 * varying layer velocity" from Fang and Chen, 2019 but a specialization for a constant layer
 * velocity. Furthermore only a first order Newton iteration is applied instead of a second order
 * using the Taylor series. This might converge slower but lead to a more stable iteration without
 * errors due to negative values in a sqrt.
 */
slowness_t TwoPointRayTracing::trace(position_t source, position_t receiver, double accuracy,
                                     int max_iterations) {
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
    b[0] = b[source_index];
    msg(b);
    msg(z);
    // eq. A5
    // TODO Many arrays only depend on the source depth. For my use case, this will be constant so
    //  they could be precalculated.
    array_t epsilon = b * b;
    msg(epsilon);
    // eq. A7
    auto h = array_t(num_layers + 1);
    for (size_t k = 1; k <= num_layers; k++) {
        h[k] = b[k] * (z[k] - z[k - 1]);
    }
    h[0] = b[source_index] * (source_z - z[source_index - 1]);
    msg(h);
    bool source_below_receiver = source_z > receiver_z;
    auto mu = array_t(num_layers + 1);
    for (size_t k = 0; k <= num_layers; k++) {
        mu[k] = mu_k(k, source_index, num_layers, source_below_receiver);
    }
    msg(mu);

    double vM = highest_velocity_between(source_z, receiver_z, model);
    msg(vM);
    // eq. A11
    array_t h_tilde = mu * h / vM;
    msg(h_tilde);

    // eq. A12
    array_t epsilon_tilde = 1 - epsilon / (vM * vM);
    msg(epsilon_tilde);

    // eq. C13
    array_t d1_unsummed = h_tilde;
    msg(d1_unsummed);
    auto d1 = nansum(d1_unsummed);
    msg(d1);
    // eq. C14
    // TODO fix delta/guards here
    array_t c0_unsummed = h_tilde / std::sqrt(epsilon_tilde);
    // this line simulates multiplication by delta epsilon but deals with the inf that would happen
    c0_unsummed[epsilon_tilde == 0] = 0;
    msg(c0_unsummed);
    auto c0 = nansum(c0_unsummed);
    msg(c0);

    // eq. C16
    auto cminus1 = 0.;
    msg(cminus1);

    // eq. C17
    array_t cminus2_unsummed = h_tilde * -0.5 / std::pow(epsilon_tilde, 1.5);
    cminus2_unsummed[epsilon_tilde == 0] = 0;
    msg(cminus2_unsummed);
    auto cminus2 = nansum(cminus2_unsummed);
    msg(cminus2);

    auto X = horizontal_distance(source_x, source_y, receiver_x, receiver_y);
    msg(X);

    // eq. C15
    auto c1 = nansum(h_tilde[epsilon_tilde == 0]);
    msg(c1);

    auto q = initial_q(X, d1, c0, c1, cminus1, cminus2);
    msg(q);

    double q_next;
    msg("Before loop");
    for (auto i = 0; i < max_iterations; ++i) {
        q_next = next_q(q, X, epsilon_tilde, h_tilde);
        msg(q_next);
        if (not std::isfinite(q_next)) {
            break;
        }
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