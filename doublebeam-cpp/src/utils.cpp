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

    std::vector<int> ray_code_to_layer_indices(const std::string& ray_code, double pz_initial,
                                               int start_index, bool include_start) {
        std::vector<int> indices;
        if (include_start) {
            indices.push_back(start_index);
        }
        bool ray_down = seismo::ray_direction_down(pz_initial);
        for (auto c : ray_code) {
            if (c == 'T') {
                start_index += ray_down ? 1 : -1;
            } else {
                ray_down = !ray_down;
            }
            indices.push_back(start_index);
        }
        return indices;
    }

    int next_layer_index(int current_index, double pz, char wave_type) {
        if (wave_type == 'R') {
            return current_index;
        } else if (wave_type == 'T') {
            return current_index + (seismo::ray_direction_down(pz) ? 1 : -1);
        } else {
            throw std::invalid_argument("Wave type not recognized: " + std::string(1, wave_type));
        }
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

    xt::xtensor<double, 2> generate_vector_arc(int num_values, double central_direction_x,
                                               double central_direction_y) {
        xt::xtensor<double, 2> angles = xt::linspace<double>(0, math::radians(180), num_values);
        auto angle_against_xaxis =
            math::angle_clockwise(1., 0., central_direction_x, central_direction_y);
        // rotate 90° to the left
        angles = angles + angle_against_xaxis - math::radians(90);
        auto x_values = xt::cos(angles);
        auto y_values = -xt::sin(angles);
        auto vectors = xt::empty<double>({num_values, 2});
        xt::view(vectors, xt::all(), xt::keep(0)) = x_values;
        xt::view(vectors, xt::all(), xt::keep(1)) = y_values;
        return vectors;
    }
}