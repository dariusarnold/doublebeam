#ifndef DOUBLEBEAM_CPP_TIMING_HPP
#define DOUBLEBEAM_CPP_TIMING_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>


struct TimingResults {
    // all time results are in nanoseconds
    double mean;
    double standard_deviation;
    uint64_t number_of_runs;
};

std::ostream& operator<<(std::ostream& os, const TimingResults& results);

template <typename InputIterator>
std::pair<typename InputIterator::value_type, typename InputIterator::value_type>
mean_and_standard_deviation(InputIterator first, InputIterator last) {
    double mean = std::accumulate(first, last, 0.) / std::distance(first, last);
    double sum = 0;
    std::for_each(first, last, [&](double x) { sum += (x - mean) * (x - mean); });
    return {mean, std::sqrt(sum / (std::distance(first, last) - 1))};
}

template <uint64_t RunTimeMilliSeconds = 4000, typename F, typename... Args>
TimingResults measure_runtime(F func, Args&&... args) {
    auto maxRunTime = std::chrono::milliseconds(RunTimeMilliSeconds);
    std::vector<double> runtimes;
    std::chrono::system_clock::time_point b;
    auto start_time = std::chrono::high_resolution_clock::now();
    do {
        auto a = std::chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        b = std::chrono::high_resolution_clock::now();
        runtimes.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count());
    } while (std::chrono::duration_cast<std::chrono::milliseconds>(b - start_time) <= maxRunTime);
    auto [mean, std_deviation] = mean_and_standard_deviation(runtimes.begin(), runtimes.end());
    return {mean, std_deviation, runtimes.size()};
}

#endif // DOUBLEBEAM_CPP_TIMING_HPP
