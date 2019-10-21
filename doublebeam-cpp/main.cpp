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
    auto source_beam_centres = grid_coordinates(400, 500, 400, 500, 0, 2, 2);
    auto targets = grid_coordinates(400, 500, 400, 500, 450, 2, 2);
    FractureParameters fractures(400, 1, 0, 61, 40, 120, 41);
    db.algorithm(source_beam_centres, targets,
                 SeismoData("/home/darius/masterarbeit/output_0degrees"), fractures, 10, 40, 0.006);
    return 0;
}
