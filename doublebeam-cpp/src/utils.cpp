#include <algorithm>
#include <xtensor/xview.hpp>

#include "utils.hpp"


namespace seismo {

    bool ray_direction_down(double pz) {
        return pz > 0;
    }


    std::tuple<double, double, double> slowness_3D(double theta, double phi, double velocity) {
        auto px = 1. / velocity * sin(theta) * cos(phi);
        auto py = 1. / velocity * sin(theta) * sin(phi);
        auto pz = 1. / velocity * cos(theta);
        return {px, py, pz};
    }

    std::vector<int> ray_code_to_layer_indices(const std::vector<WaveType>& ray_code,
                                               double pz_initial, int start_index,
                                               bool include_start) {
        std::vector<int> indices;
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

    std::vector<int> ray_code_to_layer_indices(const std::string& code, double pz_initial,
                                               int start_index, bool include_start) {
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
            if (c == 'T') {
                ray_code.push_back(WaveType::Transmitted);
            } else if (c == 'R') {
                ray_code.push_back(WaveType::Reflected);
            } else {
                throw std::invalid_argument("Invalid char " + std::string(1, c) + " in ray code.");
            }
        }
        return ray_code;
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
        std::generate(result.begin(), result.end(), [&]() { return start + i++ * distance;});
        return result;
    }
} // namespace math

namespace impl {

    impl::Formatter::operator std::string() const {
        return stream.str();
    }

    CommaSeparated::operator std::string() const {
        auto s = stream.str();
        return s.substr(0, s.size() - 2);    }
} // namespace impl
