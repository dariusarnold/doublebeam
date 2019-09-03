#ifndef DOUBLEBEAM_CPP_TIMING_H
#define DOUBLEBEAM_CPP_TIMING_H


#include <cstdint>
#include <ostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>
#include <numeric>
#include <fstream>


namespace chrono = std::chrono;


struct TimingResults {
    // all time results are in nanoseconds
    double mean;
    double standard_deviation;
    uint64_t number_of_runs;
};


std::ostream& operator<<(std::ostream& os, const TimingResults& results);


template<typename InputIterator>
std::pair<typename InputIterator::value_type, typename InputIterator::value_type>
calculate_mean_and_standard_deviation(InputIterator first, InputIterator last) {
    double mean = std::accumulate(first, last, 0.) / std::distance(first, last);
    double sum = 0;
    std::for_each(first, last, [&](double x) { sum += (x - mean) * (x - mean); });
    return {mean, std::sqrt(sum / (std::distance(first, last) - 1))};
}


template<uint64_t RunTimeMilliSeconds = 4000, typename F, typename... Args>
TimingResults measure_runtime(F func, Args&& ... args) {
    auto maxRunTime = std::chrono::milliseconds(RunTimeMilliSeconds);
    std::vector<double> runtimes;
    chrono::system_clock::time_point b;
    auto start_time = chrono::high_resolution_clock::now();
    do {
        auto a = chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        b = chrono::high_resolution_clock::now();
        runtimes.push_back(chrono::duration_cast<chrono::nanoseconds>(b - a).count());
    } while (chrono::duration_cast<chrono::milliseconds>(b - start_time) <= maxRunTime);
    auto[mean, std_deviation] = calculate_mean_and_standard_deviation(runtimes.begin(), runtimes.end());
    return {mean, std_deviation, runtimes.size()};
}

#endif //DOUBLEBEAM_CPP_TIMING_H
