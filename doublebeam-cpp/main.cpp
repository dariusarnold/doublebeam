#include "doublebeam.hpp"
#include "fft.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "timing/timing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"
#include <chrono>
#include <complex>
#include <io.hpp>


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
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto db = DoubleBeam(vm);
    auto source_beam_centres = seismo::grid_coordinates(400, 500, 400, 500, 0, 1, 1);
    auto targets = seismo::grid_coordinates(400, 500, 400, 500, 450, 1, 1);
    FractureParameters fractures(400, 1, 0, 2, 40, 120, 3);
    auto result = db.algorithm(source_beam_centres, targets,
                               SeismoData("/home/darius/masterarbeit/output_0degrees"), fractures,
                               10, 40, 0.006);
    std::ofstream file{"result.txt"};
    if (file.is_open()) {
        file << result.data;
    }
    return 0;
}
