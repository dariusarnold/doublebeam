#ifndef DOUBLEBEAM_CPP_IO_HPP
#define DOUBLEBEAM_CPP_IO_HPP

#include <vector>
#include <filesystem>

#include "raytracing_types.hpp"

/**
 * Return receiver positions as an array.
 * @param path Path to receiver file.
 * @return
 */
std::vector<position_t> read_receiverfile(std::filesystem::path path);

#endif // DOUBLEBEAM_CPP_IO_HPP
