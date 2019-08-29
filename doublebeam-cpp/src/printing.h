#ifndef DOUBLEBEAM_CPP_PRINTING_H
#define DOUBLEBEAM_CPP_PRINTING_H

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> a) {
    os << "(";
    for (auto& i: a) {
        os << i << " ";
    }
    os << ")";
    return os;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    for (auto& i: v) {
        os << i << "\n";
    }
    return os;
}


template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
    os << p.first << " " << p.second;
    return os;
}


#endif //DOUBLEBEAM_CPP_PRINTING_H
