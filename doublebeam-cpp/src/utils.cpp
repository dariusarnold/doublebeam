#include <algorithm>

#include "utils.hpp"


namespace seismo {

    std::tuple<double, double, double> slowness_3D(double theta, double phi, double velocity) {
        auto px = 1. / velocity * sin(theta) * cos(phi);
        auto py = 1. / velocity * sin(theta) * sin(phi);
        auto pz = 1. / velocity * cos(theta);
        return {px, py, pz};
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

    std::vector<position_t> grid_coordinates(double x0, double x1, double y0, double y1,
                                             double depth, size_t num_x, size_t num_y) {
        double x_stepsize = num_x == 1 ? x1 - x0 : (x1 - x0) / (num_x - 1);
        double y_stepsize = num_y == 1 ? y1 - y0 : (y1 - y0) / (num_y - 1);
        std::vector<position_t> points;
        for (auto ix = 0UL; ix < num_x; ++ix) {
            for (auto iy = 0UL; iy < num_y; ++iy) {
                points.emplace_back(x0 + ix * x_stepsize, y0 + iy * y_stepsize, depth);
            }
        }
        return points;
    }
} // namespace seismo

namespace math {
    double angle(double x1, double y1, double z1, double x2, double y2, double z2, bool acute) {
        double angle = std::acos(std::clamp(
            dot(x1, y1, z1, x2, y2, z2) / (length(x1, y1, z1) * length(x2, y2, z2)), -1., 1.));
        if (acute) {
            return angle;
        }
        return M_PI - angle;
    }

    double length(double x, double y, double z) {
        return std::sqrt(x * x + y * y + z * z);
    }

    double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
        return x1 * x2 + y1 * y2 + z1 * z2;
    }

    std::vector<Vector2> generate_vector_arc(int num_values, double central_direction_x,
                                             double central_direction_y) {
        auto angle_against_xaxis =
            math::angle_clockwise(1., 0., central_direction_x, central_direction_y);
        auto angles = math::linspace(angle_against_xaxis - radians(90),
                                     angle_against_xaxis + math::radians(90), num_values);
        std::vector<Vector2> v(num_values);
        std::transform(angles.begin(), angles.end(), v.begin(), [](double angle) {
            return Vector2{std::cos(angle), -std::sin(angle)};
        });
        return v;
    }


    std::vector<double> linspace(double start, double stop, size_t num) {
        if (num == 0) {
            return {};
        }
        if (num == 1) {
            return {start};
        }
        auto distance = (stop - start) / (num - 1);
        std::vector<double> result(num);
        size_t i = 0;
        std::generate(result.begin(), result.end(), [&]() { return start + i++ * distance; });
        return result;
    }

    std::tuple<double, double, double> cross(double x1, double y1, double z1, double x2, double y2,
                                             double z2) {
        auto x = y1 * z2 - z1 * y2;
        auto y = z1 * x2 - y1 * z2;
        auto z = x1 * y2 - y1 * x2;
        return {x, y, z};
    }

    std::tuple<double, double, double> normalize(double x, double y, double z) {
        auto l = math::length(x, y, z);
        return {x / l, y / l, z / l};
    }

    bool between(double a, double x, double b) {
        if (a < b) {
            return a <= x and x <= b;
        } else{
            return b <= x and x <= b;
        }
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
