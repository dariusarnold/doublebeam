#include "io.hpp"
#include <algorithm>
#include <filesystem>
#include <iostream>


using SeismogramFileContent = std::pair<std::vector<double>, std::vector<double>>;

void save_as_binary(std::filesystem::path path, const SeismogramFileContent& seismogram) {
    save_binary(seismogram.first, path.replace_extension("t"));
    save_binary(seismogram.second, path.replace_extension("x"));
}


void convert_all_to_binary(const std::filesystem::path& project_dir) {
    if (not std::filesystem::exists(project_dir)) {
        throw std::invalid_argument("Could not find shotdata folder " + project_dir.string() + ".");
    }
    if (not std::filesystem::is_directory(project_dir)) {
        throw std::invalid_argument(project_dir.string() + " is not a directory.");
    }
    std::cout << "Converting all files in " << project_dir << "\n";
    for (auto dir_entry : std::filesystem::recursive_directory_iterator(project_dir)) {
        if (dir_entry.is_directory()) {
            continue;
        }
        auto s = read_seismogram(dir_entry);
        save_as_binary(dir_entry, s);
    }
}


int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cerr << "Specify path to shotdata folder as command line argument.";
        exit(-1);
    }
    if (argc > 2) {
        std::cerr << "Too many arguments.";
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
