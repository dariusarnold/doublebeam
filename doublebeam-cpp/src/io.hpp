#ifndef DOUBLEBEAM_CPP_IO_HPP
#define DOUBLEBEAM_CPP_IO_HPP

#include <filesystem>
#include <fstream>
#include <vector>

#include <gsl/span>

#include "seismodata.hpp"


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

#endif // DOUBLEBEAM_CPP_IO_HPP
