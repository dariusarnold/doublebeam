#ifndef DOUBLEBEAM_CPP_IO_HPP
#define DOUBLEBEAM_CPP_IO_HPP

#include <filesystem>
#include <fstream>
#include <vector>

#include <gsl/span>

#include "raytracing_types.hpp"
#include "seismodata.hpp"
#include "utils.hpp"


/**
 * Return receiver positions as an array.
 * @param path Path to receiver file.
 * @return
 */
std::vector<Receiver> read_receiverfile(std::filesystem::path path);

/**
 * Return source positions as an array.
 * @param path
 * @return
 */
std::vector<Source> read_sourcefile(std::filesystem::path path);


/**
 * Save vector in binary format.
 * @tparam T
 * @param vec
 * @param path
 */
template <typename T>
void save_binary(const std::vector<T>& vec, std::filesystem::path path) {
    std::ofstream file{path, std::ios::out | std::ios::binary};
    file.write(reinterpret_cast<const char*>(&vec[0]), vec.size() * sizeof(T));
}


/**
 * Save all seismograms in one binary file.
 * For every seismogram save timesteps, amplitudes.
 * @param seismograms Vector of seismograms.
 * @param path Filepath under which the seismograms should be stored.
 */
void save_binary_seismograms(
    const std::vector<std::pair<std::vector<double>, std::vector<double>>>& seismograms,
    const std::filesystem::path& path);

/**
 *
 * @param N Number of seismograms in file.
 * @param path
 * @return
 */
std::vector<std::pair<std::vector<double>, std::vector<double>>>
load_binary_seismograms(size_t N, const std::filesystem::path& path);


/**
 * Loads only amplitude data, not timesteps from binary file.
 * @param number_of_seismograms Number of seismograms in file.
 * @param path
 * @return
 */
void load_binary_seismograms2(const std::filesystem::path& path, size_t number_of_seismograms,
                              gsl::span<double> amplitudes);


/**
 * Load vector from binary format. Be sure to specify the correct type as a template paramter.
 * @tparam T Type with which the data was saved.
 * @param path
 * @return
 */
template <typename T>
std::vector<T> load_binary(std::filesystem::path path) {
    std::ifstream file{path, std::ios::in | std::ios::binary};
    auto N = std::filesystem::file_size(path);
    std::vector<T> vec(N / sizeof(double));
    file.read(reinterpret_cast<char*>(&vec[0]), vec.size() * sizeof(double));
    return vec;
}


/**
 * Read only amplitude data from seismogram.
 * @param path
 * @return
 */
std::vector<double> read_amplitude(std::filesystem::path path);

/**
 * Read only timesteps from seismogram
 * @param path
 * @return
 */
std::vector<double> read_timesteps(std::filesystem::path path);


/**
 *
 * @param path
 * @return Timesteps as first element, amplitude as second element.
 */
std::pair<std::vector<double>, std::vector<double>> read_seismogram(std::filesystem::path path);

/**
 * Extract number from string after equal sign.
 * @tparam Number int, double
 * @param string String with format a = N where a is any identifier (name) and N is a number.
 * @return Number that comes after the equal sign.
 */
template <typename Number>
Number extract_number(const std::string& string) {
    // +1 to exclude = from substring
    auto equal_position = string.find('=') + 1;
    char* end;
    if constexpr (std::is_floating_point_v<Number>) {
        Number res = std::strtod(string.data() + equal_position, &end);
        return res;
    }
    if constexpr (std::is_integral_v<Number>) {
        const int base = 10;
        Number res = std::strtol(string.data() + equal_position, &end, base);
        return res;
    }
}

const auto extract_double = extract_number<double>;
const auto extract_int = extract_number<int>;

inline std::filesystem::path extract_path(const std::string& string) {
    // +1 to exclude = from substring
    auto equal_position = string.find('=') + 1;
    // could be made more efficient, doesnt matter for parsing the file.
    auto part_after_equal = string.substr(equal_position);
    auto first_not_whitespace = part_after_equal.find_first_not_of(' ');
    auto path_string = part_after_equal.substr(first_not_whitespace);
    std::filesystem::path p(path_string);
    return p;
}

inline bool contains(const std::string& string, std::string_view search_term) {
    auto n = string.find(search_term);
    return n != std::string::npos;
}

struct SourceBeamCenterParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        while (std::getline(config_file, line)) {
            if (contains(line, "x0")) {
                x0.get() = extract_double(line);
                continue;
            }
            if (contains(line, "x1")) {
                x1.get() = extract_double(line);
                continue;
            }
            if (contains(line, "y0")) {
                y0.get() = extract_double(line);
                continue;
            }
            if (contains(line, "y1")) {
                y1.get() = extract_double(line);
                continue;
            }
            if (contains(line, "z")) {
                z.get() = extract_double(line);
                continue;
            }
            if (contains(line, "num_x")) {
                num_x = extract_int(line);
                continue;
            }
            if (contains(line, "num_y")) {
                num_y = extract_int(line);
                continue;
            }
            break;
        }
    }
    Meter x0, x1, y0, y1, z;
    int num_x, num_y;

    //    SourceBeamCenterParams(std::ifstream& config_file) {}
};

struct FractureParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        while (std::getline(config_file, line)) {
            if (contains(line, "phi_hat_x")) {
                phi_hat.x = extract_double(line);
                continue;
            }
            if (contains(line, "phi_hat_y")) {
                phi_hat.y = extract_double(line);
                continue;
            }
            if (contains(line, "num_orientations")) {
                num_orientations = extract_int(line);
                continue;
            }
            if (contains(line, "spacing_min")) {
                spacings_min = extract_double(line);
                continue;
            }
            if (contains(line, "spacing_max")) {
                spacings_max = extract_double(line);
                continue;
            }
            if (contains(line, "num_spacings")) {
                num_spacings = extract_int(line);
                continue;
            }
            break;
        }
    }
    math::Vector2 phi_hat;
    int num_orientations;
    double spacings_min, spacings_max;
    int num_spacings;
};

struct BeamParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        while (std::getline(config_file, line)) {
            if (contains(line, "width")) {
                width = extract_double(line);
                continue;
            }
            if (contains(line, "frequency")) {
                frequency = extract_double(line);
                continue;
            }
            if (contains(line, "window_length")) {
                window_length = extract_double(line);
                continue;
            }
            if (contains(line, "max_stacking_distance")) {
                max_stacking_distance = extract_double(line);
                continue;
            }
            break;
        }
    }
    double width;
    double frequency;
    double window_length;
    double max_stacking_distance;
};

struct VelocityModelParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        while (std::getline(config_file, line)) {
            if (contains(line, "file")) {
                path = extract_path(line);
                break;
            }
        }
    }
    std::filesystem::path path;
};

struct SeismoDataParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        while (std::getline(config_file, line)) {
            if (contains(line, "path")) {
                path = extract_path(line);
                continue;
            }
            if (contains(line, "source_frequency")) {
                source_frequency = Frequency(extract_double(line));
                continue;
            }
            break;
        }
    }
    std::filesystem::path path;
    Frequency source_frequency;
};

struct TargetParams {
    inline void init_from_file(std::ifstream& config_file) {
        std::string line;
        double x = 0, y = 0, z = 0;
        while (std::getline(config_file, line)) {
            if (contains(line, "x")) {
                x = extract_number<double>(line);
                continue;
            }
            if (contains(line, "y")) {
                y = extract_number<double>(line);
                continue;
            }
            if (contains(line, "z")) {
                z = extract_number<double>(line);
                continue;
            }
            break;
        }
        position = {Meter(x), Meter(y), Meter(z)};
    }
    Position position;
};

struct DoubleBeamOptions {
    VelocityModelParams model_params{};
    SeismoDataParams seismo_data_params{};
    TargetParams target{};
    SourceBeamCenterParams sbc_params{};
    FractureParams fracture_params{};
    BeamParams beam_params{};
};

inline DoubleBeamOptions read_config_file(const std::filesystem::path& config_path) {
    std::ifstream config_file{config_path};
    if (not config_file.is_open()) {
        throw std::runtime_error("Failed to open config file " + config_path.string());
    }
    std::string line;
    DoubleBeamOptions options;
    while (std::getline(config_file, line)) {
        if (contains(line, "[model]")) {
            options.model_params.init_from_file(config_file);
        }
        if (contains(line, "[data]")) {
            options.seismo_data_params.init_from_file(config_file);
        }
        if (contains(line, "[target]")) {
            options.target.init_from_file(config_file);
        }
        if (contains(line, "[source beam centers]")) {
            options.sbc_params.init_from_file(config_file);
        }
        if (contains(line, "[fractures]")) {
            options.fracture_params.init_from_file(config_file);
        }
        if (contains(line, "[beam]")) {
            options.beam_params.init_from_file(config_file);
        }
    }
    return options;
}


#endif // DOUBLEBEAM_CPP_IO_HPP
