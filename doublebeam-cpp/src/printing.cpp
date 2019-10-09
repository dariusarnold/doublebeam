#include "printing.hpp"


std::ostream& operator<<(std::ostream& os, const WaveType& wave_type) {
    return os << to_char(wave_type);
}

std::ostream& operator<<(std::ostream& os, std::vector<WaveType> v) {
    for (auto wave_type : v) {
        os << wave_type;
    }
    return os;
}
