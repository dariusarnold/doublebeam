#include <fstream>
#include <optional>
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

std::vector<Receiver> read_receiverfile(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file) {
        throw std::runtime_error(impl::Formatter() << "Couldn't open receiver  file " << path);
    }
    auto n = 0;
    // read number of receiver positions in file
    file >> n;
    size_t index;
    std::vector<Receiver> receivers;
    receivers.reserve(n);
    double x, y, z;
    for (auto i = n; i > 0; --i) {
        file >> index >> x >> y >> z;
        receivers.push_back({x, y, z, index});
    }
    return receivers;
}


std::vector<Source> read_sourcefile(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file) {
        throw std::runtime_error(impl::Formatter() << "Couldn't open source file " << path);
    }
    std::string int_regex{R"((\d+))"};
    std::string float_regex{R"((\d+(?:\.\d*)?|\.\d+))"};
    std::regex nsrc_regex{R"(\s*nsrc\s*=\s*(\d+))"};
    std::regex coordinates_regex{R"(\s*xsource\s*=\s*)" + float_regex + R"([\n\s]+ysource\s*=\s*)" +
                                 float_regex + R"([\n\s]+zsource\s*=\s*)" + float_regex};
    std::regex source_regex{R"(source\s*=\s*)" + int_regex};
    // find number of sources in file
    std::string line;
    std::smatch nsources_match;
    std::optional<int> nsrc = -1;
    while (std::getline(file, line)) {
        if (std::regex_search(line, nsources_match, nsrc_regex)) {
            nsrc = std::stoi(nsources_match[1]);
            break;
        }
    }
    if (not nsrc) {
        // did not find nrsc in file (sentinel value is the same)
        throw std::runtime_error(
            impl::Formatter() << "Number of sources (nsrc = ...) not specified in file " << path);
    }
    std::vector<Source> sources;
    sources.reserve(nsrc.value());
    double x, y, z;
    // read lines until a full match (x, y and z) is found
    std::string multilines;
    size_t index = 1;
    std::smatch index_match;
    while (std::getline(file, line)) {
        multilines += line + "\n";
        if (std::regex_search(multilines, nsources_match, coordinates_regex) and
            std::regex_search(multilines, index_match, source_regex)) {
            x = std::stod(nsources_match[1]);
            y = std::stod(nsources_match[2]);
            z = std::stod(nsources_match[3]);
            index = std::stoll(index_match[1]);
            sources.push_back({x, y, z, index});
            multilines.clear();
        }
    }
    if (sources.size() != nsrc) {
        throw std::runtime_error(impl::Formatter()
                                 << "Source file " << path << " malformed. Specified "
                                 << nsrc.value() << " sources, contains only " << sources.size());
    }
    return sources;
}

std::vector<double> read_column(std::filesystem::path path, int column) {
    if (column < 0 or column > 1) {
        throw std::runtime_error(impl::Formatter() << "Invalid column index " << column
                                                   << ", allowed values are 0, 1.");
    }
    std::ifstream file{path};
    if (not file.is_open()) {
        throw std::runtime_error(impl::Formatter()
                                 << "Couldn't open file " << path << " for reading data.");
    }
    constexpr const char comment = '#';
    std::string s;
    std::vector<double> vec;
    std::array<double, 2> p;
    while (std::getline(file, s)) {
        // read lines until first line not starting with comment character
        if (s[s.find_first_not_of(' ')] != comment) {
            // push back this line
            std::istringstream(s) >> p[0] >> p[1];
            vec.push_back(p[column]);
            break;
        }
    }
    // then push back other lines
    while (file >> p[0] >> p[1]) {
        vec.push_back(p[column]);
    }
    return vec;
}

std::vector<double> read_amplitude(std::filesystem::path path) {
    return read_column(path, 1);
}

std::vector<double> read_timesteps(std::filesystem::path path) {
    return read_column(path, 0);
}
