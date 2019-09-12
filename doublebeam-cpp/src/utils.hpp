#ifndef DOUBLEBEAM_CPP_UTILS_HPP
#define DOUBLEBEAM_CPP_UTILS_HPP

#include <cmath>
#include <tuple>
#include <type_traits>


namespace seismo {

    /**
     * Calculate slowness from velocity and angles.
     * For geometric definitions see chapter 3.2.1 in Cerveny - Seismic ray theory
     * @param theta Angle against downgoing vertical axis (z) in rad, increasing
     * upwards. 0 <= theta <= pi.
     * @param phi Angle against horizontal x-axis in rad, increasing towards y-axis.
     * 0 <= phi <= 2pi.
     * @param velocity Velocity in m/s.
     * @return Slowness px, py, pz
     */
    std::tuple<double, double, double> slowness_3D(double theta, double phi, double velocity);


    /**
     * Return true if ray is going downwards, false otherwise.
     * @param pz Vertical component of slowness vector.
     */
    bool ray_direction_down(double pz);
} // namespace seismo

namespace math {

    /**
     * Convert radians to degrees.
     * @tparam T
     * @param degrees Angle in degrees.
     * @return Angle in radians.
     */
    template <typename T>
    double radians(T degrees) {
        constexpr double factor = M_PI / 180.;
        if constexpr (std::is_integral<T>::value) {
            degrees = static_cast<double>(degrees);
        }
        return degrees * factor;
    }


    /**
     * Check if two numbers have the same sign.
     * Zero is treated as a positive number.
     * @tparam T Number type.
     * @return True if the sign of both numbers is the same, false otherwise.
     */
    template <typename T>
    bool same_sign(T a, T b) {
        return std::signbit(a) == std::signbit(b);
    }

    template <typename T>
    /**
     * Calculate clockwise angle in radians from the first vector to the second vector in the x_y
     * horizontal plane.
     * @tparam T floating point type
     * @param x1 x component of first vector.
     * @param y1 y component of first vector.
     * @param x2 x component of second vector.
     * @param y2 y component of second vector.
     * @return
     */
    double angle_clockwise(T x1, T y1, T x2, T y2) {
        auto angle1 = std::atan2(x1, y1);
        auto angle2 = std::atan2(x2, y2);
        auto difference = angle2 - angle1;
        if (difference < 0)
            return 2 * M_PI + difference;
        return difference;
    }

} // namespace math
#endif // DOUBLEBEAM_CPP_UTILS_HPP
