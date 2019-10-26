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
    while (file >> index >> x >> y >> z) {
        receivers.push_back({x, y, z, index});
    }
    return receivers;
}


std::string stream_to_string(std::istream& is) {
    std::ostringstream buffer;
    buffer << is.rdbuf();
    return buffer.str();
}


std::vector<Source> read_sourcefile(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file) {
        throw std::runtime_error(impl::Formatter() << "Couldn't open source file " << path);
    }
    std::string int_regex{R"((\d+))"};
    std::string float_regex{R"((\d+(?:\.\d*)?|\.\d+))"};
    std::regex nsrc_regex{R"(\s*nsrc\s*=\s*)" + int_regex};
    std::regex source_regex{R"(\bsource\s*=\s*)" + int_regex + R"([\n\s]+xsource\s*=\s*)" +
                            float_regex + R"([\n\s]+ysource\s*=\s*)" + float_regex +
                            R"([\n\s]+zsource\s*=\s*)" + float_regex};
    // find number of sources in file
    std::string line;
    std::smatch nsources_match;
    std::optional<int> nsrc;
    while (std::getline(file, line)) {
        if (std::regex_search(line, nsources_match, nsrc_regex)) {
            nsrc = std::stoi(nsources_match[1]);
            break;
        }
    }
    if (not nsrc) {
        // did not find nrsc in file (optional value is still unset)
        throw std::runtime_error(
            impl::Formatter() << "Number of sources (nsrc = ...) not specified in file " << path);
    }
    std::vector<Source> sources;
    sources.reserve(nsrc.value());
    auto filecontent = stream_to_string(file);
    auto sources_begin = std::sregex_iterator(filecontent.begin(), filecontent.end(), source_regex);
    auto sources_end = std::sregex_iterator();
    for (auto it = sources_begin; it != sources_end; ++it) {
        auto match = *it;
        size_t index = std::stoll(match[1]);
        double x = std::stod(match[2]);
        double y = std::stod(match[3]);
        double z = std::stod(match[4]);
        sources.push_back({x, y, z, index});
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
    size_t pos{0};
    while (std::getline(file, s)) {
        p[0] = std::stod(s, &pos);
        p[1] = std::stod(s.substr(pos, s.size()));
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

std::pair<std::vector<double>, std::vector<double>> read_seismogram(std::filesystem::path path) {
    std::ifstream file{path};
    if (not file.is_open()) {
        throw std::runtime_error(impl::Formatter() << " Failed opening file " << path);
    }
    std::string s;
    constexpr char comment = '#';
    double t, x;
    std::vector<double> times, amplitdues;
    while (std::getline(file, s)) {
        // read lines until first line not starting with comment character
        if (s[s.find_first_not_of(' ')] != comment) {
            // push back this line
            std::istringstream(s) >> t >> x;
            times.push_back(t);
            amplitdues.push_back(x);
            break;
        }
    }
    while (file >> t >> x) {
        times.push_back(t);
        amplitdues.push_back(x);
    }
    return {times, amplitdues};
}

void save_binary_seismograms(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& seismograms,
    const std::filesystem::path& path) {
    std::ofstream file{path, std::ios::binary};
    for (const auto& seismogram : seismograms) {
        file.write((const char*)(seismogram.first.data()),
                   seismogram.first.size() * sizeof(double));
        file.write((const char*)(seismogram.second.data()),
                   seismogram.second.size() * sizeof(double));
    }
}

[[nodiscard]] std::vector<std::pair<std::vector<double>, std::vector<double>>>
load_binary_seismograms(size_t N, const std::filesystem::path& path) {
    std::vector<std::pair<std::vector<double>, std::vector<double>>> seismograms;
    std::ifstream file{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error("Failed to open file " + path.string());
    }
    // number of entries in the times/amplitudes vector of a single seismogram
    auto size_seismogram = std::filesystem::file_size(path) / N / 2 / sizeof(double);
    std::vector<double> t(size_seismogram), x(size_seismogram);
    for (auto i = 0U; i < N; ++i) {
        file.read(reinterpret_cast<char*>(t.data()), size_seismogram * sizeof(double));
        file.read(reinterpret_cast<char*>(x.data()), size_seismogram * sizeof(double));
        seismograms.emplace_back(t, x);
    }
    return seismograms;
}