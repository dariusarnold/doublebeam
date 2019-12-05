#include <chrono>
#include <complex>

#include "doublebeam.hpp"
#include "eigen_helpers.hpp"
#include "io.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
#include "units.hpp"
#include "utils.hpp"


std::ostream& operator<<(std::ostream& os, const RaySegment& segment) {
    os << "Raysegment: " << segment.data.front() << " to " << segment.data.back() << " ( "
       << segment.data.size() << ") steps";
    return os;
}


std::ostream& operator<<(std::ostream& os, const Ray& ray) {
    auto i = 0;
    for (auto& segment : ray.segments) {
        std::cout << i << " " << segment << "\n";
        ++i;
    }
    return os;
}


int main() {
    std::ios_base::sync_with_stdio(false);
    auto vm = read_velocity_file("/home/darius/git/doublebeam/doublebeam-cpp/velocity_model.txt");
    auto db = DoubleBeam(vm);
    auto source_beam_centres = seismo::grid_coordinates(4500, 8500, 4500, 8500, 0, 4, 4);
    auto targets = seismo::grid_coordinates(6500, 6500, 6500, 6500, 2400, 1, 1);
    // TODO fracture depth is unused
    FractureParameters fractures(2400, 0, 1, 61, 100, 300, 41);
    auto data = SeismoData("/home/darius/masterarbeit/output_0degrees");
    auto a = std::chrono::high_resolution_clock::now();
    auto result = db.algorithm(source_beam_centres, targets, data, fractures, 244_meter,
                               40_angular_from_hertz, 0.06);
    auto b = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime db : " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count()
              << " s" << std::endl;
    std::string result_filename = "result";
    std::filesystem::path result_path(result_filename + ".txt");
    int index = 0;
    while (std::filesystem::exists(result_path)) {
        ++index;
        result_path = std::filesystem::path(result_filename + (std::to_string(index)) + ".txt");
    }
    std::ofstream file{result_path};
    if (file.is_open()) {
        file << result.data;
    }
    return 0;
}