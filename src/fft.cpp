#include <iostream>

#include <gsl/span>

#include "fft.hpp"


PlanCache::PlanCache(size_t initial_cache_size) : plans() {
    plans.reserve(initial_cache_size);
}

Plan::Plan(std::vector<double>& in) :
        N(in.size()),
        out(N / 2 + 1),
        p(fftw_plan_dft_r2c_1d(N, in.data(), reinterpret_cast<fftw_complex*>(out.data()),
                               FFTW_MEASURE)) {}

Plan::Plan(gsl::span<const double> in) :
        N(in.size()),
        out(N / 2 + 1),
        p(fftw_plan_dft_r2c_1d(N, const_cast<double*>(in.data()),
                               reinterpret_cast<fftw_complex*>(out.data()), FFTW_MEASURE)) {}

std::vector<std::complex<double>, Allocator<std::complex<double>>>& Plan::execute() {
    fftw_execute(p.get());
    return out;
}

FFT::cvector FFT::execute(std::vector<double>& in) {
    if (in.empty()) {
        return {};
    }
    auto& plan = plans.get_plan(in);
    const auto& res = plan.execute();
    return std::vector(res.begin(), res.end());
}

std::vector<FFT::cdouble, Allocator<std::complex<double>>>
FFT::execute(gsl::span<const double> in) {
    if (in.empty()) {
        return {};
    }
    auto& plan = plans.get_plan(in);
    return plan.execute();
}
