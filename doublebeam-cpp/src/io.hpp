#ifndef DOUBLEBEAM_CPP_IO_HPP
#define DOUBLEBEAM_CPP_IO_HPP

#include <filesystem>
#include <fstream>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>
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
 * Loads only amplitude data, not timesteps from binary file.
 * @param number_of_seismograms Number of seismograms in file.
 * @param path
 * @return
 */
void load_binary_seismograms(const std::filesystem::path& path, size_t number_of_seismograms,
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


struct SourceBeamCenterParams {
    Meter x0, x1, y0, y1, z;
    int num_x, num_y;

    friend std::ostream& operator<<(std::ostream& os, const SourceBeamCenterParams& sbc_params) {
        fmt::print(os, "[source_beam_centers]\n");
        fmt::print(os, "x0 = {}\nx1 = {}\ny0 = {}\ny1 = {}\nz = {}\nnum_x = {}\nnum_y = {}\n\n",
                   sbc_params.x0.get(), sbc_params.x1.get(), sbc_params.y0.get(),
                   sbc_params.y1.get(), sbc_params.z.get(), sbc_params.num_x, sbc_params.num_y);
        return os;
    }
};

struct FractureParams {
    math::Vector2 phi_hat{1, 0};
    int num_orientations;
    Meter spacings_min, spacings_max;
    int num_spacings;

    friend std::ostream& operator<<(std::ostream& os, const FractureParams& params) {
        fmt::print(os, "[fractures]\n");
        fmt::print(os,
                   "phi_hat_x = {}\nphi_hat_y = {}\nnum_orientations = {}\nspacing_min = "
                   "{}\nspacing_max = {}\nnum_spacings = {}\n\n",
                   params.phi_hat.x, params.phi_hat.y, params.num_orientations,
                   params.spacings_min.get(), params.spacings_max.get(), params.num_spacings);
        return os;
    }
};

struct BeamParams {
    Meter width;
    Frequency reference_frequency;
    Frequency source_frequency;
    Second window_length;
    Meter max_stacking_distance;

    friend std::ostream& operator<<(std::ostream& os, const BeamParams& beam_params) {
        fmt::print(os, "[beam]\n");
        fmt::print(os,
                   "width = {}\nreference_frequency = {}\nwindow_length = "
                   "{}\nmax_stacking_distance = {}\nsource_frequency={}\n\n",
                   beam_params.width.get(), beam_params.reference_frequency.get(),
                   beam_params.window_length.get(), beam_params.max_stacking_distance.get(),
                   beam_params.source_frequency.get());
        return os;
    }
};

struct SeismoDataParams {
    std::filesystem::path data_path;
    std::filesystem::path velocity_model_path;

    friend std::ostream& operator<<(std::ostream& os, const SeismoDataParams& data_params) {
        fmt::print(os, "[data]\n");
        fmt::print(os, "path = {}\n", data_params.data_path.string());
        fmt::print(os, "model = {}\n\n", data_params.velocity_model_path.string());
        return os;
    }
};

struct TargetParams {
    Position position;

    friend std::ostream& operator<<(std::ostream& os, const TargetParams& target_params) {
        fmt::print(os, "[target]\n");
        fmt::print(os, "x = {}\ny = {}\nz = {}\n\n", target_params.position.x.get(),
                   target_params.position.y.get(), target_params.position.z.get());
        return os;
    }
};


struct DoubleBeamOptions {
    SeismoDataParams seismo_data_params{};
    TargetParams target{};
    SourceBeamCenterParams sbc_params{};
    FractureParams fracture_params{};
    BeamParams beam_params{};

    inline friend std::ostream& operator<<(std::ostream& os, const DoubleBeamOptions& options) {
        os << options.seismo_data_params << options.target << options.sbc_params
           << options.fracture_params << options.beam_params;
        return os;
    }
};

#endif // DOUBLEBEAM_CPP_IO_HPP
