#ifndef DOUBLEBEAM_CPP_IO_HPP
#define DOUBLEBEAM_CPP_IO_HPP

#include <filesystem>
#include <fstream>
#include <vector>

#include "raytracing_types.hpp"

/**
 * Return receiver positions as an array.
 * @param path Path to receiver file.
 * @return
 */
std::vector<position_t> read_receiverfile(std::filesystem::path path);

/**
 * Return source positions as an array.
 * @param path
 * @return
 */
std::vector<position_t> read_sourcefile(std::filesystem::path path);


/**
 * Save vector in binary format.
 * @tparam T
 * @param vec
 * @param path
 */
template <typename T>
void save_binary(const std::vector<T>& vec, std::filesystem::path path) {
    std::ofstream file{path, std::ios::out | std::ios::binary};
    file.write((const char*)(&vec[0]), vec.size() * sizeof(T));
}

/**
 * Load vector from binary format. Be sure to specify the correct type as a template paramter.
 * @tparam T Type with which the data was saved.
 * @param path
 * @return
 */
template <typename T>
std::vector<T> load_binary(std::filesystem::path path) {
    std::ifstream file{path, std::ios::in | std::ios::binary};
    file.seekg(0, file.end);
    auto N = file.tellg();
    file.seekg(0, file.beg);
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

#endif // DOUBLEBEAM_CPP_IO_HPP
