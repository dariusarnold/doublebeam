#include <fstream>
#include <regex>

#include "io.hpp"
#include "utils.hpp"

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
        throw std::runtime_error(impl::Formatter() << "Couldn't open receiver  file " << path);
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


std::vector<position_t> read_sourcefile(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file) {
        throw std::runtime_error(impl::Formatter() << "Couldn't open source file " << path);
    }
    std::string int_regex{R"((\d+))"};
    std::regex nsrc_regex{R"(nsrc\s+=\s+(\d+))"};
    std::string float_regex{R"((\d+(?:\.\d*)?|\.\d+))"};
    std::regex source_regex{R"(xsource\s+=\s+)" + float_regex + R"([\n\s]+ysource\s+=\s+)" +
                            float_regex + R"([\n\s]+zsource\s+=\s+)" + float_regex};
    // find number of sources in file
    std::string line;
    std::smatch nsources_match;
    auto nsrc = 0UL;
    while (std::getline(file, line)) {
        if (std::regex_search(line, nsources_match, nsrc_regex)) {
            nsrc = std::stoi(nsources_match[1]);
            break;
        }
    }
    std::vector<position_t> source_positions;
    source_positions.reserve(nsrc);
    double x, y, z;
    // read lines until a full match (x, y and z) is found
    std::string multilines;
    while (std::getline(file, line)) {
        multilines += line + "\n";
        if (std::regex_search(multilines, nsources_match, source_regex)) {
            x = std::stod(nsources_match[1]);
            y = std::stod(nsources_match[2]);
            z = std::stod(nsources_match[3]);
            source_positions.emplace_back(x, y, z);
            multilines.clear();
        }
    }
    if (source_positions.size() != nsrc) {
        throw std::runtime_error(impl::Formatter()
                                 << "Source file " << path << " malformed. Specified " << nsrc
                                 << " sources, contains only " << source_positions.size());
    }
    return source_positions;
}
