#ifndef DOUBLEBEAM_CPP_UTILS_HPP
#define DOUBLEBEAM_CPP_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <vector>

#include "raytracing_types.hpp"


// Helper code used for implementation, not for user code
namespace impl {

    // Use to throw exceptions with customized error message by streaming into a Formatter instance.
    class Formatter {
    public:
        Formatter(const std::string& sep = "");

        template <typename T>
        Formatter& operator<<(const T& t) {
            if (not first_used) {
                // On first use, only stream t to avoid prepending a separator
                first_used = true;
                stream << t;
            } else {
                // stream separator and then new value so there is never a trailing separator
                stream << separator << t;
            }
            return *this;
        }

        /**
         * Alllow streaming Formatter to stream by converting it explicitly to a string.
         */
        friend std::ostream& operator<<(std::ostream& os, const Formatter& f);

        operator std::string() const;

    private:
        bool first_used = false;
        std::stringstream stream{};
        std::string separator;
    };

    // Base case: T is an arithmetic type (a floating point type), so just use T.
    template <bool B, typename T>
    struct get_type_helper {
        using type = T;
    };

    // Specialized type for the case where T is not an arithmetic type.
    // This is covers std::complex<T>, which is a class.
    template <typename T>
    struct get_type_helper<false, T> {
        using type = typename T::value_type;
    };

    // The structs typedef named type will be of type T if T is an arithmetic type. Otherwise it
    // will be of type T::value_type. The second case handles std::complex numbers.
    // Motivation:
    // Goertzel algorithm works with complex and real input data. In both cases, a single complex
    // value is returned. This struct extracts the std::complex's value type if T is complex or just
    // uses T as a type if T is a basic floating point type. I had to use this to declare the return
    // type of the templated Goertzel algorithm, which should be a complex type of the base floating
    // point type. Not using this would result in a std::complex<std::complex<double>> when passing
    // in complex data.
    template <typename T>
    struct value_type_or_type {
        using type = typename get_type_helper<std::is_arithmetic_v<T>, T>::type;
    };

} // namespace impl


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
    constexpr bool ray_direction_down(double pz) {
        return pz > 0;
    }

    /**
     * Generate ray code from string.
     * Raise invalid_argument exception when string is not valid (contains characters other than T
     * and R).
     * @param s String of T, R, where T specifies the transmitted wave and R the reflected wave.
     * @return Ray code.
     */
    std::vector<WaveType> make_ray_code(const std::string s);

    /**
     * Transform a ray code to a sequence of layer indices.
     * @param ray_code Specifies which ray type (reflected, transmitted) to take at an interface.
     * Valid values are 'R' for reflected and 'T' for transmitted.
     * @param pz_initial Initial slowness of the ray, used to get initial direction. A negative
     * value represents a upgoing ray, a positive value represents a downgoing ray.
     * @param start_index At which layer index the sequence should start.
     * @return
     */
    // TODO this function doesn't work for turning rays.
    std::vector<std::ptrdiff_t> ray_code_to_layer_indices(const std::vector<WaveType>& ray_code,
                                                          double pz_initial,
                                                          std::ptrdiff_t start_index = 0,
                                                          bool include_start = true);

    std::vector<std::ptrdiff_t> ray_code_to_layer_indices(const std::string& code,
                                                          double pz_initial,
                                                          std::ptrdiff_t start_index = 0,
                                                          bool include_start = true);

    /**
     * Calculate next layer index for a given layer index and state at an interface (state before
     * crossing the interface).
     * @param current_index Index of current layer.
     * @param pz Vertical slowness at the interface, before crossing or reflection.
     * @param wave_type
     * @return
     */
    int next_layer_index(int current_index, double pz, WaveType wave_type);

    /**
     * Generate coordinates of a evenly spaced grid.
     * @param x0 X start coordinate.
     * @param x1 X end coordinate.
     * @param y0 Y start coordinate.
     * @param y1 Y end coordinate.
     * @param depth Depth of all grid points.
     * @param num_x Number of points along the x axis.
     * @param num_y Number of points along the y axis.
     * @return Container of all grid locations with shape (num_x, num_y, 3) where the first axis
     * is the x axis, the second axis is the y axis and the third axis contains the z values.
     */
    std::vector<position_t> grid_coordinates(double x0, double x1, double y0, double y1,
                                             double depth, size_t num_x, size_t num_y);
} // namespace seismo


namespace math {

    /**
     * Convert radians to degrees.
     * @tparam T
     * @param degrees Angle in degrees.
     * @return Angle in radians.
     */
    template <typename T>
    constexpr double radians(T degrees) {
        constexpr double factor = M_PI / 180.;
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
     * @param x1 X component of first vector.
     * @param y1 X component of first vector.
     * @param z1 Z component of first vector.
     * @param x2 X component of second vector.
     * @param y2 Y component of second vector.
     * @param z2 Z component of second vector.
     */
    double dot(double x1, double y1, double z1, double x2, double y2, double z2);

    /**
     * Cross product of two 3D vectors.
     * @param x1 X component of first vector.
     * @param y1 X component of first vector.
     * @param z1 Z component of first vector.
     * @param x2 X component of second vector.
     * @param y2 Y component of second vector.
     * @param z2 Z component of second vector.
     * @return
     */
    std::tuple<double, double, double> cross(double x1, double y1, double z1, double x2, double y2,
                                             double z2);


    /**
     * Scale vector to new length while keeping the direction.
     * @param vector x,y, z tuple
     * @param new_length New length of vector.
     * @return Vector with length new_length and direction of vector.
     */
    std::tuple<double, double, double>
    scale_vector(const std::tuple<double, double, double>& vector, double new_length);

    /**
     * Normalize vector by making it's length 1 while keeping the direction.
     * @param x
     * @param y
     * @param z
     * @return
     */
    std::tuple<double, double, double> normalize(double x, double y, double z);

    /**
     * Return determinant of 3x3 matrix
     *       |a b c|
     * |A| = |d e f|
     *       |g h i|
     */
    template <typename T>
    T det(T a, T b, T c, T d, T e, T f, T g, T h, T i) {
        return a * e * i + b * f * g + c * d * h - g * e * c - h * f * a - i * d * b;
    }

    /**
     * Return inverse of 3x3 matrix
     *         a b c ^-1
     * A^-1 = (d e f)
     *         g h i
     */
    template <typename T>
    std::tuple<T, T, T, T, T, T, T, T, T> inv(T a, T b, T c, T d, T e, T f, T g, T h, T i) {
        auto dd = 1 / det(a, b, c, d, e, f, g, h, i);
        return {dd * (e * i - f * h), dd * (c * h - b * i), dd * (b * f - c * e),
                dd * (f * g - d * i), dd * (a * i - c * g), dd * (c * d - a * f),
                dd * (d * h - e * g), dd * (b * g - a * h), dd * (a * e - b * d)};
    }

    /**
     * Right multiply matrix by column vector.
       // TODO I need a linear algebra library.
     * @tparam T
     * @param matrix
     * @param vector
     * @return
     */
    template <typename T>
    std::tuple<T, T, T> dot(std::tuple<T, T, T, T, T, T, T, T, T> matrix,
                            std::tuple<T, T, T> vector) {
        auto [a, b, c, d, e, f, g, h, i] = matrix;
        auto [x, y, z] = vector;
        return {a * x + b * y + c * z, d * x + e * y + f * z, g * x + h * y + i * z};
    }

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
        std::vector<content_t> result;
        result.reserve(length);
        if (not std::isnan(initial)) {
            result.push_back(initial);
        }
        content_t running_sum = 0.;
        for (auto a = ybegin, b = ybegin + 1; b != yend; ++a, ++b) {
            running_sum += distance * 0.5 * (*a + *b);
            result.push_back(running_sum);
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
        std::vector<content_t> result;
        result.reserve(length);
        if (not std::isnan(initial)) {
            result.push_back(initial);
        }
        content_t running_sum = 0.;
        for (auto ya = ybegin, yb = ybegin + 1, xa = xbegin, xb = xbegin + 1; yb != yend;
             ++ya, ++yb, ++xa, ++xb) {
            running_sum += (*xb - *xa) * 0.5 * (*ya + *yb);
            result.push_back(running_sum);
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

    struct Vector2 {
        double x, y;
    };

    /**
     * Generate num_values of unit vectors in a 180° arc around the central direction.
     * The first vector of the returned range will have a clockwise angle of -90° of the central
     * direction, the last vector will have a 90° clockwise angle to the central direction.
     * @param num_values The first dimension of the returned result will have this size.
     * @param central_direction_x X component of central direction.
     * @param central_direction_y Y component of central direction.
     * @return
     */
    std::vector<Vector2> generate_vector_arc(int num_values, double central_direction_x,
                                             double central_direction_y);

    /**
     * Create vector filled with num values from start to stop.
     * @param start
     * @param stop
     * @param num
     * @return
     */
    std::vector<double> linspace(double start, double stop, size_t num = 50);


    /**
     * Calculate value of a single frequency bin for the discrete Fourier transform.
     * @tparam T Scalar data type. Can be complex or real.
     * @param data Input sequence of length N.
     * @param target_frequency_bin Frequencies are restricted to the frequency bins of the discrete
     * Fourier transform (DFT). The frequency to be evaluated is given as target_frequency = 2 * pi
     * * k/N k is an index number going from 0...N-1
     * @return Value of DFT frequency bin.
     */
    template <typename T>
    std::complex<typename impl::value_type_or_type<T>::type>
    goertzel(const std::vector<T>& data, unsigned long long target_frequency_bin) {
        if (data.size() == 0) {
            throw std::invalid_argument("Data for Goertzel algorithm is empty.");
        }
        if (target_frequency_bin > data.size())
            throw std::invalid_argument(impl::Formatter()
                                        << "Invalid input bin " << target_frequency_bin << ".");
        // This will be a float type
        using type = typename impl::value_type_or_type<T>::type;
        auto N_float = static_cast<type>(data.size());
        const auto pi = static_cast<type>(M_PIl);
        const type target_frequency = 2. * pi * (target_frequency_bin / N_float);
        T s_prev_prev = data[0];
        T s_prev = data[1] + 2 * std::cos(target_frequency) * s_prev_prev;
        T s = 0;
        for (auto i = 2U; i < data.size(); ++i) {
            s = data[i] + static_cast<type>(2.) * std::cos(target_frequency) * s_prev - s_prev_prev;
            s_prev_prev = s_prev;
            s_prev = s;
        }
        return std::exp(std::complex<type>(0, 2) * pi * (target_frequency_bin / N_float)) * s_prev -
               s_prev_prev;
    }

    /**
     * Returns value of closest frequency bin to a given frequency of the discrete Fourier
     * transform.
     * @tparam T
     * @param data
     * @param frequency Target frequency, closest frequency bin of the DFT result to this frequency
     * is returned.
     * @param sampling_period Sampling period (time between two samples) of data in seconds.
     * @return
     */
    template <typename T>
    std::complex<typename impl::value_type_or_type<T>::type>
    fft_closest_frequency(const std::vector<T>& data, T frequency, double sampling_period) {
        auto sampling_frequency_Hz = 1 / sampling_period;
        auto bin = std::llround(data.size() * frequency / sampling_frequency_Hz);
        return goertzel(data, bin);
    }

    /**
     * Return true if x is between a and b (inclusive).
     */
    bool between(double a, double x, double b);

} // namespace math
#endif // DOUBLEBEAM_CPP_UTILS_HPP
