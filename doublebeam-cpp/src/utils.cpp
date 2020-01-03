#include <algorithm>

#include "utils.hpp"


namespace seismo {

    Slowness slowness_3D(Radian theta, Radian phi, Velocity velocity) {
        InverseVelocity px(1. / velocity.get() * sin(theta.get()) * cos(phi.get()));
        InverseVelocity py(1. / velocity.get() * sin(theta.get()) * sin(phi.get()));
        InverseVelocity pz(1. / velocity.get() * cos(theta.get()));
        return Slowness(px, py, pz);
    }

    std::vector<std::ptrdiff_t> ray_code_to_layer_indices(const std::vector<WaveType>& ray_code,
                                                          double pz_initial,
                                                          std::ptrdiff_t start_index,
                                                          bool include_start) {
        std::vector<std::ptrdiff_t> indices;
        if (include_start) {
            indices.push_back(start_index);
        }
        bool ray_down = seismo::ray_direction_down(pz_initial);
        for (auto c : ray_code) {
            if (c == WaveType::Transmitted) {
                start_index += ray_down ? 1 : -1;
            } else {
                ray_down = !ray_down;
            }
            indices.push_back(start_index);
        }
        return indices;
    }

    std::vector<std::ptrdiff_t> ray_code_to_layer_indices(const std::string& code,
                                                          double pz_initial,
                                                          std::ptrdiff_t start_index,
                                                          bool include_start) {
        return ray_code_to_layer_indices(make_ray_code(code), pz_initial, start_index,
                                         include_start);
    }

    int next_layer_index(int current_index, double pz, WaveType wave_type) {
        if (wave_type == WaveType::Reflected) {
            return current_index;
        } else if (wave_type == WaveType::Transmitted) {
            return current_index + (seismo::ray_direction_down(pz) ? 1 : -1);
        } else {
            throw std::invalid_argument("Wave type not recognized.");
        }
    }

    std::vector<WaveType> make_ray_code(const std::string s) {
        std::vector<WaveType> ray_code;
        for (char c : s) {
            ray_code.push_back(to_wavetype(c));
        }
        return ray_code;
    }

    std::vector<Position> grid_coordinates(Meter x0, Meter x1, Meter y0, Meter y1, Meter depth,
                                           size_t num_x, size_t num_y) {
        Meter x_stepsize = num_x == 1 ? x1 - x0 : (x1 - x0) / (num_x - 1);
        Meter y_stepsize = num_y == 1 ? y1 - y0 : (y1 - y0) / (num_y - 1);
        std::vector<Position> points;
        for (auto ix = 0UL; ix < num_x; ++ix) {
            for (auto iy = 0UL; iy < num_y; ++iy) {
                points.emplace_back(x0 + ix * x_stepsize, y0 + iy * y_stepsize, depth);
            }
        }
        return points;
    }
} // namespace seismo

namespace math {
    Radian angle(double x1, double y1, double z1, double x2, double y2, double z2, bool acute) {
        Radian angle(std::acos(std::clamp(
            dot(x1, y1, z1, x2, y2, z2) / (length(x1, y1, z1) * length(x2, y2, z2)), -1., 1.)));
        if (acute) {
            return angle;
        }
        return Radian(M_PI) - angle;
    }

    double length(double x, double y, double z) {
        return std::sqrt(x * x + y * y + z * z);
    }

    double length(const Slowness& slowness) {
        return length(slowness.px.get(), slowness.py.get(), slowness.pz.get());
    }

    double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
        return x1 * x2 + y1 * y2 + z1 * z2;
    }

    std::vector<Vector2> generate_vector_arc(int num_values, math::Vector2 central_direction) {
        auto angle_against_xaxis =
            math::angle_clockwise(1., 0., central_direction.x, central_direction.y);
        auto angles = math::linspace(angle_against_xaxis - radians(90_deg).get(),
                                     angle_against_xaxis + radians(90_deg).get(), num_values);
        std::vector<Vector2> v(num_values);
        std::transform(angles.begin(), angles.end(), v.begin(), [](double angle) {
            return Vector2{std::cos(angle), -std::sin(angle)};
        });
        return v;
    }

    std::tuple<double, double, double> cross(double x1, double y1, double z1, double x2, double y2,
                                             double z2) {
        auto x = y1 * z2 - z1 * y2;
        auto y = z1 * x2 - x1 * z2;
        auto z = x1 * y2 - y1 * x2;
        return {x, y, z};
    }

    std::tuple<double, double, double> normalize(double x, double y, double z) {
        auto l = math::length(x, y, z);
        return {x / l, y / l, z / l};
    }

    std::tuple<double, double, double>
    scale_vector(const std::tuple<double, double, double>& vector, double new_length) {
        auto [x, y, z] = vector;
        double length = math::length(x, y, z);
        return {x * new_length / length, y * new_length / length, z * new_length / length};
    }

    size_t calc_frequency_bin(size_t N, AngularFrequency target_frequency,
                              AngularFrequency sampling_frequency) {
        // bin index is frequency divided by bin width rounded
        // bin width = sampling_frequency / N (number of time domain data points)
        return std::llround(N * target_frequency.get() / sampling_frequency.get());
    }

} // namespace math

namespace impl {

    impl::Formatter::operator std::string() const {
        return stream.str();
    }

    Formatter::Formatter(const std::string& sep) : separator(sep) {}

    std::ostream& operator<<(std::ostream& os, const Formatter& f) {
        return os << std::string(f);
    }
} // namespace impl
