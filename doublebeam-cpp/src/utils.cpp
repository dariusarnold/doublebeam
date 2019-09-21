#include <algorithm>

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
                                               int start_index) {
        std::vector<int> indices{start_index};
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
}