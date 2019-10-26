#ifndef DOUBLEBEAM_CPP_FFT_HPP
#define DOUBLEBEAM_CPP_FFT_HPP

#include <algorithm>
#include <complex>
#include <cstdint>
#include <limits>
#include <new>
#include <vector>

#include <fftw3.h>


/**
 * Custom allocator using fftw_malloc and fftw_free to comply with
 * fftws alignment requirements. This enables us to use std::vector
 * of c++ complex type.
 * @tparam T
 */
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
        return reinterpret_cast<T*>(fftw_malloc(n * sizeof(T)));
    }
    void deallocate(T* p, std::size_t) noexcept {
        fftw_free(p);
    }
};


struct Plan {
    // move constructor is required because the fftw_plan frees the memory associated with the
    // input/output arrays when it is destroyed. The move constructor will invalidate the other
    // plan, so the arrays don't get double freed when the vector of plans reallocates.
    Plan(Plan&& other);
    Plan(std::vector<double>& in);
    ~Plan();
    std::vector<std::complex<double>, Allocator<std::complex<double>>>& execute();

    size_t N;
    std::vector<std::complex<double>, Allocator<std::complex<double>>> out;
    fftw_plan p;
};


class PlanCache {
public:
    PlanCache(size_t initial_cache_size);
    /**
     * Return a plan from cache.
     * if the requested plan size is not already in the cache, add it.
     * @param in
     * @return
     */
    Plan& get_plan(std::vector<double>& in);

private:
    // since I dont expect to have many different plans, a vector should be enough and a map is
    // not required, even though the vector has only a linear search to find the correct plan size.
    std::vector<Plan> plans;
};


class FFT {
public:
    FFT(size_t initial_cache_size = 2) : plans(initial_cache_size) {}
    using cdouble = std::complex<double>;
    using cvector = std::vector<cdouble>;
    cvector execute(std::vector<double>& in);
private:
    PlanCache plans;
};

#endif // DOUBLEBEAM_CPP_FFT_HPP
