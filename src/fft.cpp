/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
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
