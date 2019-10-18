#include "fft.hpp"
#include <iostream>

Plan& PlanCache::get_plan(std::vector<double>& in) {
    auto N = in.size();
    auto plan =
        std::find_if(plans.begin(), plans.end(), [N](const Plan& plan) { return plan.N == N; });
    if (plan != plans.end()) {
        return *plan;
    }
    Plan p{in};
    plans.emplace_back(std::move(p));
    return plans.back();
}

Plan::Plan(std::vector<double>& in) : N(in.size()), out(N / 2 + 1) {
    p = fftw_plan_dft_r2c_1d(N, in.data(), reinterpret_cast<fftw_complex*>(out.data()),
                             FFTW_MEASURE);
}

Plan::~Plan() {
    fftw_destroy_plan(p);
}

std::vector<std::complex<double>, Allocator<std::complex<double>>>& Plan::execute() {
    fftw_execute(p);
    return out;
}

FFT::cvector FFT::execute(std::vector<double>& in) {
    auto plan = plans.get_plan(in);
    return plan.execute();
}
