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
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <regex>

#include "io.hpp"


namespace fs = std::filesystem;

using SeismogramFileContent = std::pair<std::vector<double>, std::vector<double>>;


std::vector<fs::path> get_all_directories(const fs::path& path) {
    std::vector<fs::path> paths;
    std::copy_if(fs::directory_iterator(path), fs::directory_iterator(), std::back_inserter(paths),
                 [](auto& dir_entry) { return dir_entry.is_directory(); });
    return paths;
}


std::vector<fs::path> get_all_files(const fs::path& path, const std::regex& seismogram_file) {
    std::vector<fs::path> files;
    auto filter = [&](auto& dir_entry) {
        return (dir_entry.is_regular_file() and
                std::regex_search(dir_entry.path().filename().string(), seismogram_file));
    };
    std::copy_if(fs::directory_iterator(path), fs::directory_iterator(), std::back_inserter(files),
                 filter);
    return files;
}


void convert_content_to_binary(const fs::path& source_dir, const std::regex& seismogram_file) {
    std::vector<SeismogramFileContent> seismograms;
    auto files = get_all_files(source_dir, seismogram_file);
    if (files.empty()) {
        throw std::length_error("Found no seismogram files in " + source_dir.string());
    }
    std::sort(files.begin(), files.end());
    std::cout << "Converting " << files.front().string() << " to " << files.back().string() << "\n";
    std::transform(files.begin(), files.end(), std::back_inserter(seismograms), read_seismogram);
    std::cout << "Saved " << (source_dir / config::get_binary_seismogram_filename()).string()
              << "\n";
    save_binary_seismograms(seismograms, source_dir / config::get_binary_seismogram_filename());
}


void convert_all_to_binary(const std::filesystem::path& shotdata_path,
                           const std::regex& seismogram_file) {
    if (not std::filesystem::exists(shotdata_path)) {
        throw std::invalid_argument("Could not find shotdata folder " + shotdata_path.string() +
                                    ".");
    }
    if (not std::filesystem::is_directory(shotdata_path)) {
        throw std::invalid_argument(shotdata_path.string() + " is not a directory.");
    }
    std::cout << "Converting all files in " << shotdata_path << std::endl;
    auto paths = get_all_directories(shotdata_path);
    std::sort(paths.begin(), paths.end());
    auto i = 0;
    for (const auto& source_dir : paths) {
        std::cout << ++i << "/" << paths.size() << ": " << source_dir << std::endl;
        convert_content_to_binary(source_dir, seismogram_file);
    }
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Specify path to shotdata folder as first positional argument and seismogram "
                     "file regex as second positional argument.\n";
        exit(-1);
    }
    std::ios::sync_with_stdio(false);
    std::filesystem::path p(argv[1]);
    std::regex seismogram_file(argv[2]);
    auto a = std::chrono::high_resolution_clock::now();
    convert_all_to_binary(p, seismogram_file);
    auto b = std::chrono::high_resolution_clock::now();
    std::cout << "Took " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count()
              << " s\n";
}
