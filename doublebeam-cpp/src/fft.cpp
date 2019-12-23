#include "fft.hpp"
#include <iostream>

Plan& PlanCache::get_plan(std::vector<double>& in) {
    auto N = in.size();
    auto plan = std::find_if(plans.begin(), plans.end(), [N](const Plan& p) { return p.N == N; });
    if (plan != plans.end()) {
        return *plan;
    }
    return plans.emplace_back(in);
}

PlanCache::PlanCache(size_t initial_cache_size) : plans() {
    plans.reserve(initial_cache_size);
}

Plan::Plan(std::vector<double>& in) :
        N(in.size()),
        out(N / 2 + 1),
        p(fftw_plan_dft_r2c_1d(N, in.data(), reinterpret_cast<fftw_complex*>(out.data()),
                               FFTW_MEASURE)) {}


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
