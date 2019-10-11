#include "doublebeam.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "timing/timing.hpp"
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

int nmain() {
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto tp = TwoPointRayTracing(vm);
    for (int i = 0; i < 10; ++i) {
        auto a = std::chrono::high_resolution_clock::now();
        tp.trace({50, 10, i*2}, {1, 10+i, 450});
        auto b = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count()
                  << std::endl;
    }
    return 0;
}


int main() {
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto db = DoubleBeam(vm);
    auto sources = grid_coordinates(400, 500, 400, 500, 0, 2, 2);
    auto targets = grid_coordinates(400, 500, 400, 500, 450, 2, 2);
    FractureParameters fractures(400, 1, 0, 61, 40, 120, 41);
    db.algorithm(sources, targets, fractures, 10, 40, 0.006);
    return 0;
}