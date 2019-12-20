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

#define USEDEBUG false
#if USEDEBUG
#define msg(x) std::cout << #x << " = " << x << std::endl
#else
#define msg(x)
#endif

/**
 * Sum values in array treating nan as zero.
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


std::string stringify(position_t pos) {
    auto [x, y, z] = pos;
    return "(" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ")";
}

// eq. A3
template <typename A>
TwoPointRayTracing::array_t X_tilde(double q, const A& epsilon_tilde, const A& h_tilde) {
    return h_tilde / std::sqrt(1 / (q * q) + epsilon_tilde);
}

// eq. B9
template <typename A>
TwoPointRayTracing::array_t X_tilde_prime(double q, const A& epsilon_tilde, const A& h_tilde) {
    return h_tilde / std::pow(1 + epsilon_tilde * q * q, 1.5);
}

// eq. B10
template <typename A>
TwoPointRayTracing::array_t X_tilde_double_prime(double q, const A& epsilon_tilde,
                                                 const A& h_tilde) {
    return -3 * h_tilde * epsilon_tilde * q / std::pow(1 + epsilon_tilde * q * q, 2.5);
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

template <typename A>
double f_tilde_double_prime(double q, const A& epsilon_tilde, const A& h_tilde) {
    return nansum(X_tilde_double_prime(q, epsilon_tilde, h_tilde));
}

// Solve nonlinear equation (5) using Newton algorithm with cubic iteration.
template <typename A>
double next_q(double q, double X, const A& epsilon_tilde, const A& h_tilde) {
    double f_double_prime = f_tilde_double_prime(q, epsilon_tilde, h_tilde);
    double f_prime = f_tilde_prime(q, epsilon_tilde, h_tilde);
    double f = f_tilde(q, X, epsilon_tilde, h_tilde);
    // use approximation which is valid for small a.
    // 1-sqrt(1-a) = a/2 + a^2/8 + O(a^3)
    // see http://numbers.computation.free.fr/Constants/Algorithms/newton.html
    auto q_new = q - f / f_prime * (1 + (f * f_double_prime / (2 * f_prime * f_prime)));
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


// Calculate starting value for q using either C1 or C6 depending on value of c1.
// For constant velocity layers c_minus1 will always be 0, so it is removed from the formulas.
double initial_q(double X, double d1, double c0, double c1, double cminus2) {
    auto alpha1 = d1;
    msg(alpha1);
    auto alpha2 = safe_divide(c0 * c0 * c0, -c0 * cminus2);
    msg(alpha2);
    auto beta1 = safe_divide(d1 * cminus2, c0 * cminus2);
    msg(beta1);
    auto beta2 = safe_divide(c0 * c0, -c0 * cminus2);
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
        auto q = X / alpha1 + (X / c1 - X / alpha1 - c0 / c1) *
                                  Heaviside(X - safe_divide(alpha1 * c0, alpha1 - c1));
        return q;
    }
}

double horizontal_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
}


TwoPointRayTracing::TwoPointRayTracing(VelocityModel velocity_model) :
        model(std::move(velocity_model)) {}


Slowness TwoPointRayTracing::trace(Position source, Position receiver, double accuracy,
                                   int max_iterations) {
    auto [min_depth, max_depth] = std::minmax(source.z.get(), receiver.z.get());
    if (not model.in_model(source.x, source.y, source.z)) {
        throw std::domain_error(impl::Formatter()
                                << "Source at " << source << " outside of model.");
    }
    if (not model.in_model(receiver.x, receiver.y, receiver.z)) {
        throw std::domain_error(impl::Formatter()
                                << "Receiver at " << receiver << "  outside of model.");
    }
    auto source_index = model.layer_index(source.z).value();
    auto receiver_index = model.layer_index(receiver.z).value();
    auto [min_index, max_index] = std::minmax(source_index, receiver_index);
    auto num_layers = model.number_of_interfaces_between(source.z, receiver.z) + 1;
    msg(source_index);
    msg(receiver_index);
    msg(num_layers);
    // Only the layers between source and receiver and the layers containing source and receiver
    // matter. All others are ignored. Furthermore the depth in z will be set to source/receiver
    // depth at the top and bottom.
    array_t z(num_layers + 1);
    array_t b(num_layers);
    z[0] = min_depth;
    z[num_layers] = max_depth;
    b[0] = model[min_index].velocity.get();
    for (size_t i = 1; i < num_layers; ++i) {
        z[i] = model[min_index + i - 1].bot_depth.get();
        b[i] = model[min_index + i].velocity.get();
    }
    msg(z);
    msg(b);
    array_t epsilon = b * b;
    msg(epsilon);
    array_t h(num_layers);
    if (num_layers == 1) {
        h = b * std::abs((source.z - receiver.z).get());
    } else {
        for (size_t i = 0; i < num_layers; ++i) {
            h[i] = b[i] * (z[i + 1] - z[i]);
        }
    }
    msg(h);
    double vM = b.max();
    msg(vM);
    array_t h_tilde = h / vM;
    msg(h_tilde);
    array_t epsilon_tilde = 1 - epsilon / (vM * vM);
    msg(epsilon_tilde);
    double d1 = h_tilde.sum();
    // TODO find better way for sum, nansum evaluates into temporary array
    double c0 = nansum((h_tilde / sqrt(epsilon_tilde))[epsilon_tilde != 0]);
    double c1 = nansum(h_tilde[epsilon_tilde == 0]);
    double cn2 = -0.5 * nansum((h_tilde / pow(epsilon_tilde, 1.5))[epsilon_tilde != 0]);
    msg(d1);
    msg(c0);
    msg(c1);
    msg(cn2);

    auto X =
        horizontal_distance(source.x.get(), source.y.get(), receiver.x.get(), receiver.y.get());
    auto q = initial_q(X, d1, c0, c1, cn2);
    for (; max_iterations > 0; --max_iterations) {
        double q_next = next_q(q, X, epsilon_tilde, h_tilde);
        if (std::abs(q - q_next) < accuracy) {
            q = q_next;
            break;
        }
        q = q_next;
    }
    double horizontal_slowness = q_to_p(q, vM);
    Velocity c = model.eval_at(source.x, source.y, source.z).value();
    InverseVelocity vertical_slowness(
        std::sqrt(std::pow(c.get(), -2) - horizontal_slowness * horizontal_slowness));
    bool source_below_receiver = source.z > receiver.z;
    if (source_below_receiver) {
        // ray should travel upward from source in this case
        vertical_slowness *= -1;
    }
    // calculate angle to x axis
    double phi =
        math::angle_clockwise((receiver.x - source.x).get(), (receiver.y - source.y).get(), 1., 0.);
    // horizontal slowness is combination of px and py, now we have to split it up into its parts
    InverseVelocity px(std::cos(phi) * horizontal_slowness);
    InverseVelocity py(std::sin(phi) * horizontal_slowness);
    return {px, py, vertical_slowness};
}
