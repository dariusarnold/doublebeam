#include "doublebeam.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "timing/timing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"
#include <chrono>
#include <complex>
#include <fftw3.h>
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

template <typename T>
struct Allocator {
    typedef T value_type;
    Allocator() = default;
    template <typename U>
    constexpr Allocator(const Allocator<U>&) noexcept {}
    T* allocate(std::size_t n) {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
            throw std::bad_alloc();
        }
        return (T*)(fftw_malloc(n * sizeof(T)));
    }
    void deallocate(T* p, std::size_t) noexcept {
        fftw_free(p);
    }
};

template <typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, Alloc>& vec) {
    for (const auto& el : vec) {
        os << el << " ";
    }
    return os;
}

struct Plan {
    ~Plan() {
        fftw_destroy_plan(p);
    }
    void plan(size_t N, double* in, std::complex<double>* out) {
        p = fftw_plan_dft_r2c_1d(N, in, reinterpret_cast<fftw_complex*>(out), FFTW_MEASURE);
    }
    void execute() {
        fftw_execute(p);
    }

private:
    fftw_plan p;
};

int main() {
    const int N = 4;
    std::vector<double> in{0, 1, 2, 3};
    std::vector<std::complex<double>, Allocator<std::complex<double>>> out2(in.size() / 2 + 1);
    Plan p;
    std::cout << in << std::endl;
    std::cout << std::endl;
    p.plan(N, in.data(), out2.data());
    p.execute();
    std::cout << out2 << std::endl;
    return 0;
}
/*
int main() {
    auto vm = read_velocity_file("/home/darius/git/doublebeam/fang2019model.txt");
    auto db = DoubleBeam(vm);
    auto sources = grid_coordinates(400, 500, 400, 500, 0, 2, 2);
    auto targets = grid_coordinates(400, 500, 400, 500, 450, 2, 2);
    FractureParameters fractures(400, 1, 0, 61, 40, 120, 41);
    db.algorithm(sources, targets, fractures, 10, 40, 0.006);
    return 0;
}
*/