#ifndef DOUBLEBEAM_CPP_PRINTING_HPP
#define DOUBLEBEAM_CPP_PRINTING_HPP

#include <cstddef>
#include <ostream>
#include <array>
#include <vector>
#include <tuple>


template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> a) {
    os << "(";
    for (auto& i : a) {
        os << i << " ";
    }
    os << ")";
    return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    for (auto& i : v) {
        os << i << " ";
    }
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
    // print last element without separator
    for (auto i = std::begin(a); i != std::end(a) - 1; ++i) {
        os << *i << ", ";
    }
    os << *(std::end(a) - 1);
    return os;
}


#endif // DOUBLEBEAM_CPP_PRINTING_HPP
