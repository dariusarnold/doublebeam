#ifndef DOUBLEBEAM_CPP_UTILS_H
#define DOUBLEBEAM_CPP_UTILS_H

#include <tuple>
#include <cmath>
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

}

namespace math {

    /**
     * Convert radians to degrees.
     * @tparam T
     * @param degrees Angle in degrees.
     * @return Angle in radians.
     */
    template<typename T>
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
    template<typename T>
    bool same_sign(T a, T b) {
        return std::signbit(a) == std::signbit(b);
    }

}
#endif //DOUBLEBEAM_CPP_UTILS_H
