#include "io.hpp"
#include "utils.hpp"
#include <fstream>

std::istream& operator>>(std::istream& is, position_t& pos) {
    // skip first column (receiver number)
    int i;
    is >> i;
    // get x, y, z position
    is >> std::get<0>(pos) >> std::get<1>(pos) >> std::get<2>(pos);
    return is;
}

std::vector<position_t> read_receiverfile(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file) {
        throw std::runtime_error(impl::Formatter() << "Couldn't open file " << path);
    }
    auto n = 0;
    // read number of receiver positions in file
    file >> n;
    std::vector<position_t> receiver_positions(n);
    for (auto i = 0; i < n; ++i) {
        position_t x;
        file >> x;
        receiver_positions[i] = x;
    }
    return receiver_positions;
}
