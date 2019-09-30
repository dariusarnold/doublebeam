#ifndef DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP
#define DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP

#include <array>
#include <complex>
#include <cstddef>


using state_type = std::array<double, 7>;

using position_t = std::tuple<double, double, double>;

using slowness_t = std::tuple<double, double, double>;

using complex = std::complex<double>;

namespace Index {
    /*
     * Indices of variables in state_type.
     * X, Y, Z are cartesian coordinates.
     * PX, PY, PZ are components of slowness vector.
     * T is travel time.
     */
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t Z = 2;
    static constexpr size_t PX = 3;
    static constexpr size_t PY = 4;
    static constexpr size_t PZ = 5;
    static constexpr size_t T = 6;
}; // namespace Index

/**
 * Specifies wave type to take when crossing an interface. Can be used to build ray codes (sequences
 * of wave types), which specify the behaviour of a ray during it's travel through the velocity
 * model.
 */
enum class WaveType {
    Transmitted,
    Reflected
};

#endif // DOUBLEBEAM_CPP_RAYTRACING_TYPES_HPP
