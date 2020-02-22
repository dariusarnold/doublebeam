#ifndef DOUBLEBEAM_CPP_PRINTING_HPP
#define DOUBLEBEAM_CPP_PRINTING_HPP

#include "raytracing_types.hpp"
#include <array>
#include <cstddef>
#include <ostream>
#include <tuple>
#include <vector>


template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> a) {
    os << "(";
    for (auto& i : a) {
        os << i << " ";
    }
    os << ")";
    return os;
}

/**
 * Print wavetypes without separating spaces.
 */
std::ostream& operator<<(std::ostream& os, std::vector<WaveType> v);

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    if (v.empty()) {
        return os;
    }
    // print last element without separator
    for (auto i = v.begin(); i != v.end() - 1; ++i) {
        os << *i << ", ";
    }
    os << *(v.end() - 1);
    return os;
}


template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
    os << p.first << " " << p.second;
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T, T, T>& t) {
    auto [x, y, z] = t;
    os << x << " " << y << " " << z;
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& a) {
    if (a.size() == 0) {
        return os;
    }
    // print last element without separator
    for (auto i = std::begin(a); i != std::end(a) - 1; ++i) {
        os << *i << ", ";
    }
    os << *(std::end(a) - 1);
    return os;
}

std::ostream& operator<<(std::ostream& is, const WaveType& wave_type);

#endif // DOUBLEBEAM_CPP_PRINTING_HPP
