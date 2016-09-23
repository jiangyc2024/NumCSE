#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

// @file most_significant_bit.cpp This file contains benchmarking
// of methods to find most significant bit.

using namespace boost::accumulators;

using accumulator_t = boost::accumulators::accumulator_set<
    double,
    stats<tag::mean, tag::min, tag::max>
    >;

using time_point_t = std::chrono::high_resolution_clock::time_point;
using duration_t = std::chrono::nanoseconds;

template <class Function>
accumulator_t time(Function && F,
                   double min_runtime = 20., // minimum runtime in seconds
                   unsigned int min_iterations = 10, // repeat at least 10 times
                   unsigned int max_iterations = 10000000) { // repeat at most 1e7 times

    accumulator_t time_accumulator;

    {
        unsigned it = 0; // iteration counter

        volatile long unsigned int k = 0;
        while (it < max_iterations
               && (it < min_iterations
               || sum(time_accumulator) < min_runtime)) {

            unsigned int i = rand();

            time_point_t start = std::chrono::high_resolution_clock::now();

            k += F(i);

            duration_t elapsed = std::chrono::duration_cast<duration_t>
                (std::chrono::high_resolution_clock::now()  - start);

            time_accumulator(elapsed.count() * 1e-9);

            ++it;
        }
        std::cout << k << std::endl;
    }

    return time_accumulator;
}

template <typename T>
inline unsigned int div_find(T N) {
    T p = 2;
    unsigned int last_one = 0;
    for(unsigned int k = 0; k < sizeof(T) * 8 - 1; ++k) {
        if( N % p == 0 ) {
            last_one = k;
        }
        p = p << 1;
    }
    return last_one;
}

template <typename T>
inline unsigned int loop_find(T N) {
    T p = 1;
    unsigned int last_one = 0;
    for(unsigned int k = 0; k < sizeof(T) * 8; ++k) {
        if( (~N & p) == 0 ) {
            last_one = k;
        }
        p = p << 1;
    }
    return last_one;
}

template <typename T>
inline int float_find(T N) {
    return std::ceil(std::log2(N));
}

template <typename T>
inline int builtin_find(T N) {
    return __builtin_clzl (N);
}

int main() {

    std::vector<accumulator_t> accs;

    accs.push_back( time( div_find<long unsigned int>, 5 ) );
    accs.push_back( time( loop_find<long unsigned int>, 5 ) );
    accs.push_back( time( float_find<long unsigned int>, 5 ) );
    accs.push_back( time( builtin_find<long unsigned int>, 5 ) );

    std::cout  << std::setprecision(3) << std::scientific
               << std::setw(15) << "mean"
               << std::setw(15) << "min"
               << std::setw(15) << "max"
               << std::endl;
    for(auto &acc : accs ) {
        std::cout  << std::setprecision(3) << std::scientific
                   << std::setw(15) << mean(acc)
                   << std::setw(15) << min(acc)
                   << std::setw(15) << max(acc)
                   << std::endl;
    }
}
