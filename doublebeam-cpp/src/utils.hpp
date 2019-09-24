#ifndef DOUBLEBEAM_CPP_UTILS_HPP
#define DOUBLEBEAM_CPP_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xtensor/xtensor.hpp>


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

    /**
     * Transform a ray code to a sequence of layer indices.
     * @param ray_code Specifies which ray type (reflected, transmitted) to take at an interface.
     * Valid values are 'R' for reflected and 'T' for transmitted.
     * @param pz_initial Initial slowness of the ray, used to get initial direction. A negative
     * value represents a upgoing ray, a positive value represents a downgoing ray.
     * @param start_index At which layer index the sequence should start.
     * @return
     */
    std::vector<int> ray_code_to_layer_indices(const std::string& ray_code, double pz_initial,
                                               int start_index = 0);

    /**
     * Calculate next layer index for a given layer index and state at an interface (state before
     * crossing the interface).
     * @param current_index Index of current layer.
     * @param pz Vertical slowness at the interface, before crossing or reflection.
     * @param wave_type
     * @return
     */
    int next_layer_index(int current_index, double pz, char wave_type);
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
    template <typename T>
    double angle_clockwise(T x1, T y1, T x2, T y2) {
        auto angle1 = std::atan2(x1, y1);
        auto angle2 = std::atan2(x2, y2);
        auto difference = angle2 - angle1;
        if (difference < 0) return 2 * M_PI + difference;
        return difference;
    }

    /**
     * Calculate angle between two 3D vectors in rad.
     * @param x1
     * @param y1
     * @param z1
     * @param x2
     * @param y2
     * @param z2
     * @return
     */
    double angle(double x1, double y1, double z1, double x2, double y2, double z2,
                 bool acute = true);

    /**
     * Calculate length of 3D vector.
     * @param x
     * @param y
     * @param z
     * @return
     */
    double length(double x, double y, double z);

    /**
     * Dot product of two 3D vectors.
     */
    double dot(double x1, double y1, double z1, double x2, double y2, double z2);

    /**
     * Round value to given number of digits after the decimal dot.
     * @tparam D Floating point type.
     * @tparam I Integral type.
     * @param value Value to round.
     * @param places To how many digits the value is rounded.
     */
    template <typename D, typename I>
    D round(D value, I places) {
        return std::floor(value * std::pow(10, places) + 0.5) / std::pow(10, places);
    }


    template <typename T, typename = void>
    struct is_container : std::false_type {};

    template <typename T>
    struct is_container<T, std::conditional_t<false,
                                              std::void_t<decltype(std::declval<T>().begin()),
                                                          decltype(std::declval<T>().end())>,
                                              void>> : public std::true_type {};

    template <typename T, typename = void>
    struct is_iterator : std::false_type {};

    template <typename T>
    struct is_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category>>
            : std::true_type {};

    /**
     * Cumulatively integrate y using the composite trapezoidal rule, where distance is the equal
     * spacing between values of y.
     * @tparam Iterator
     * @param ybegin Start of sequence of y values.
     * @param yend One past the end of sequence of y values.
     * @param distance Spacing between elements of y.
     * @param initial If given, insert this value at the beginning of the returned vector. If not
     * given, result will have one less element than y.
     * @return Vector of cumulative integration results of y.
     */
    template <typename Iterator, std::enable_if_t<is_iterator<Iterator>::value, int> = 0>
    std::vector<double>
    cumtrapz(Iterator ybegin, Iterator yend, typename Iterator::value_type distance = 1.,
             typename Iterator::value_type initial =
                 std::numeric_limits<typename Iterator::value_type>::quiet_NaN()) {
        using content_t = typename Iterator::value_type;
        auto length = std::distance(ybegin, yend) - 1;
        auto i = 0;
        // if initial value was given, increase length of result by one, insert initial value and
        // start from index 1 so initial value is not overwritten. Initialize with correct size
        // before instead of inserting at zero after to avoid overhead.
        if (not std::isnan(initial)) {
            length += 1;
            i += 1;
        }
        std::vector<content_t> result(length);
        if (not std::isnan(initial)) {
            result[0] = initial;
        }
        content_t running_sum = 0.;
        for (auto a = ybegin, b = ybegin + 1; b != yend; ++a, ++b) {
            running_sum += distance * 0.5 * (*a + *b);
            result[i] = running_sum;
            ++i;
        }
        return result;
    }

    /**
     * Overload for cumtrapz that works with the container directly instead of iterators.
     * @tparam Sequence
     * @param y
     * @param distance
     * @param initial
     * @return
     */
    template <typename Sequence, std::enable_if_t<is_container<Sequence>::value, int> = 0>
    std::vector<typename Sequence::value_type>
    cumtrapz(const Sequence& y, typename Sequence::value_type distance = 1.,
             typename Sequence::value_type initial =
                 std::numeric_limits<typename Sequence::value_type>::quiet_NaN()) {
        return cumtrapz(y.begin(), y.end(), distance, initial);
    }


    /**
     * Cumulatively integrate y(x) using the composite trapezoidal rule.
     * @tparam Iterator
     * @param ybegin Start of sequence of y values.
     * @param yend One past the end of sequence of y values.
     * @param ybegin Start of sequence of x values.
     * @param yend One past the end of sequence of x values.
     * @param initial If given, insert this value at the beginning of the returned vector. If not
     * given, result will have one less element than y.
     * @return Vector of cumulative integration results of y(x).
     */
    template <typename Iterator, std::enable_if_t<is_iterator<Iterator>::value, int> = 0>
    std::vector<double>
    cumtrapz(Iterator ybegin, Iterator yend, Iterator xbegin, Iterator xend,
             typename Iterator::value_type initial =
                 std::numeric_limits<typename Iterator::value_type>::quiet_NaN()) {
        using content_t = typename Iterator::value_type;
        auto length = std::distance(ybegin, yend) - 1;
        if (std::distance(xbegin, xend) - 1 != length) {
            throw std::invalid_argument("Length of x and y not the same.");
        }
        auto i = 0;
        // if initial value was given, increase length of result by one, insert initial value and
        // start from index 1 so initial value is not overwritten. Initialize with correct size
        // before instead of inserting at zero after to avoid overhead.
        if (not std::isnan(initial)) {
            length += 1;
            i += 1;
        }
        std::vector<content_t> result(length);
        if (not std::isnan(initial)) {
            result[0] = initial;
        }
        content_t running_sum = 0.;
        for (auto ya = ybegin, yb = ybegin + 1, xa = xbegin, xb = xbegin + 1; yb != yend;
             ++ya, ++yb, ++xa, ++xb) {
            running_sum += (*xb - *xa) * 0.5 * (*ya + *yb);
            result[i] = running_sum;
            ++i;
        }
        return result;
    }

    /**
     * Overload for cumtrapz that works directly with the containers instead of iterators.
     * @tparam Sequence1
     * @tparam Sequence2
     * @param y
     * @param x
     * @param initial
     * @return
     */
    template <typename Sequence1, typename Sequence2,
              std::enable_if_t<is_container<Sequence1>::value, int> = 0,
              std::enable_if_t<is_container<Sequence2>::value, int> = 0>
    std::vector<double>
    cumtrapz(const Sequence1& y, const Sequence2& x,
             typename Sequence1::value_type initial =
                 std::numeric_limits<typename Sequence1::value_type>::quiet_NaN()) {
        return cumtrapz(y.begin(), y.end(), x.begin(), x.end(), initial);
    }


    /**
     * Generate num_values of unit vectors in a 180° arc around the central direction.
     * The first vector of the returned range will have a clockwise angle of -90° of the central
     * direction, the last vector will have a 90° clockwise angle to the central direction.
     * @param num_values The first dimension of the returned result will have this size.
     * @param central_direction_x X component of central direction.
     * @param central_direction_y Y component of central direction.
     * @return
     */
    xt::xtensor<double, 2> generate_vector_arc(int num_values, double central_direction_x,
                                               double central_direction_y);

} // namespace math
#endif // DOUBLEBEAM_CPP_UTILS_HPP
