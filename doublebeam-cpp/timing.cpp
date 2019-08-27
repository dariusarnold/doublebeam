#include <chrono>
#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <algorithm>


struct TimingResults{
    // all time units are in microseconds
    double mean;
    double standard_deviation;
    uint64_t number_of_runs;
};


template <typename T>
std::pair<T, T> calculate_mean_and_standard_deviation(std::vector<T> v){
    T mean = std::accumulate(v.begin(), v.end(), 0.) / v.size();
    T sum = 0;
    std::for_each(v.begin(), v.end(), [&](double x){sum += (x - mean) * (x - mean);});
    return {mean, std::sqrt(sum / (v.size() - 1))};
}


template<typename TimeUnit, typename F, typename... Args>
TimingResults measure_runtime(F func, Args&&... args){
    uint64_t number_of_runs = 0;
    double total_runtime_ms = 4000;
    std::vector<double> runtimes;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::chrono::system_clock::time_point b;
    do {
        auto a = std::chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        b = std::chrono::high_resolution_clock::now();
        runtimes.push_back(std::chrono::duration_cast<TimeUnit>(b-a).count());
        number_of_runs++;
    } while (std::chrono::duration_cast<std::chrono::milliseconds>(b-start_time).count() <= total_runtime_ms);
    auto [mean, std_deviation] = calculate_mean_and_standard_deviation(runtimes);
    return {mean, std_deviation, number_of_runs};
}