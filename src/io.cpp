/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include <fstream>
#include <optional>
#include <regex>

#include "io.hpp"
#include "utils.hpp"


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
    Meter x, y, z;
    while (file >> index >> x >> y >> z) {
        receivers.emplace_back(x, y, z, index);
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
        Meter x(std::stod(match[2]));
        Meter y(std::stod(match[3]));
        Meter z(std::stod(match[4]));
        sources.emplace_back(x, y, z, index);
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
        file.write(reinterpret_cast<const char*>(seismogram.first.data()),
                   seismogram.first.size() * sizeof(double));
        file.write(reinterpret_cast<const char*>(seismogram.second.data()),
                   seismogram.second.size() * sizeof(double));
    }
}


void load_binary_seismograms(const std::filesystem::path& path, size_t number_of_seismograms,
                              gsl::span<double> amplitudes) {
    std::ifstream binary_file{path, std::ios::binary};
    if (not binary_file) {
        throw std::runtime_error("Failed to open file " + path.string());
    }
    // Number of data points in the amplitude vector of a single seismograms
    const auto samples_seismogram =
        std::filesystem::file_size(path) / sizeof(double) / 2 / number_of_seismograms;
    for (auto i = 0U; i < number_of_seismograms; ++i) {
        // ignore timestep data
        binary_file.ignore(samples_seismogram * sizeof(double));
        // read amplitude data
        binary_file.read(reinterpret_cast<char*>(amplitudes.data() + i * samples_seismogram),
                         samples_seismogram * sizeof(double));
    }
}