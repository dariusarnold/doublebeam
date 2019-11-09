#include <chrono>
#include <complex>

#include "doublebeam.hpp"
#include "eigen_helpers.hpp"
#include "io.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
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
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto db = DoubleBeam(vm);
    auto source_beam_centres = seismo::grid_coordinates(450, 500, 450, 500, 0, 1, 1);
    // TODO integration doesnt stop exactly at depth, only at layers
    auto targets = seismo::grid_coordinates(420, 500, 410, 500, 450, 1, 1);
    FractureParameters fractures(400, 1, 0, 2, 40, 120, 3);
    auto a = std::chrono::high_resolution_clock::now();
    auto result = db.algorithm(source_beam_centres, targets,
                               SeismoData("/home/darius/masterarbeit/output_0degrees"), fractures,
                               10, 40, 0.006);
    auto b = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime db : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count() << " ms"
              << std::endl;
    std::ofstream file{"result.txt"};
    if (file.is_open()) {
        file << result.data;
    }
    return 0;
}