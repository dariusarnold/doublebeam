#include "io.hpp"
#include <algorithm>
#include <filesystem>
#include <iostream>


namespace fs = std::filesystem;

using SeismogramFileContent = std::pair<std::vector<double>, std::vector<double>>;


std::vector<fs::path> get_all_directories(const fs::path& path) {
    std::vector<fs::path> paths;
    std::copy_if(fs::directory_iterator(path), fs::directory_iterator(), std::back_inserter(paths),
                 [](auto& dir_entry) { return dir_entry.is_directory(); });
    return paths;
}


std::vector<fs::path> get_all_files(const fs::path& path, const std::string& file_extension = "") {
    std::vector<fs::path> files;
    auto filter = [&](auto& dir_entry) {
        return (dir_entry.is_regular_file() and
                (file_extension.empty() ? true : dir_entry.path().extension() == file_extension));
    };
    std::copy_if(fs::directory_iterator(path), fs::directory_iterator(), std::back_inserter(files),
                 filter);
    return files;
}


void convert_content_to_binary(const fs::path& source_dir) {
    std::vector<SeismogramFileContent> seismograms;
    auto files = get_all_files(source_dir, ".txt");
    if (files.empty()) {
        throw std::length_error("Found no seismogram files in " + source_dir.string());
    }
    std::sort(files.begin(), files.end());
    std::transform(files.begin(), files.end(), std::back_inserter(seismograms), read_seismogram);
    save_binary_seismograms(seismograms, source_dir / "data.bin");
}


void convert_all_to_binary(const std::filesystem::path& project_dir) {
    if (not std::filesystem::exists(project_dir)) {
        throw std::invalid_argument("Could not find shotdata folder " + project_dir.string() + ".");
    }
    if (not std::filesystem::is_directory(project_dir)) {
        throw std::invalid_argument(project_dir.string() + " is not a directory.");
    }
    std::cout << "Converting all files in " << project_dir << std::endl;
    auto paths = get_all_directories(project_dir);
    auto i = 0;
    for (const auto& source_dir : paths) {
        std::cout << ++i << "/" << paths.size() << ": " << source_dir << std::endl;
        convert_content_to_binary(source_dir);
    }
}


int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cerr << "Specify path to shotdata folder as command line argument.\n";
        exit(-1);
    }
    if (argc > 2) {
        std::cerr << "Too many arguments.\n";
        exit(-1);
    }
    std::ios::sync_with_stdio(false);
    std::filesystem::path p(argv[1]);
    auto a = std::chrono::high_resolution_clock::now();
    convert_all_to_binary(p);
    auto b = std::chrono::high_resolution_clock::now();
    std::cout << "Took " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count()
              << " s\n";
}
