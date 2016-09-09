#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

// boost accumulators
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "timer.h"
#include "strassen.cpp"

#include "polyfit.hpp"

using namespace Eigen;

// we use accumulators to collect important characteristics
//  like minimum-, maximum- and mean-runtime
// note: just collecting the mean runtime is sufficient
//  but collection of extremal values can give valuable insights
using namespace boost::accumulators;

using accumulator_t = boost::accumulators::accumulator_set<
    double, // collect number of seconds in double precision
    stats<tag::sum, tag::min, tag::max, tag::mean> // collect min, max, mean
>;

// tell the compiler that p is used in some sense and should not be
//  optimized away
template <typename T>
static void escape (T&& p) {
    asm volatile("" :  : "g"(p) : "memory");
}

VectorXd asymptotic_complexity(std::vector<unsigned> dims_, std::vector<accumulator_t>& timings, 
  unsigned offset=0) {
    assert(dims_.size() == timings.size());
    assert(dims_.size() > offset);
    // retrieve mean values
    VectorXd means(timings.size()-offset);
    for (unsigned i=0; i<timings.size()-offset; ++i) {
        means[i] = mean(timings[i+offset]);
    }
    // convert std::vector<unsigned> into Eigen::VectorXd
    VectorXd dims = Map<Matrix<unsigned, 1, -1>>(dims_.data()+offset, 1,
        dims_.size()-offset).cast<double>();
    
    return polyfit(dims.array().log().matrix(), means.array().log().matrix(), 1);
}

int main()
{
    unsigned seed = (unsigned) time(0); // store seed for rng
    
    double min_runtime = 100.; // minimum runtime in seconds
    unsigned int min_iterations = 10; // repeat at least 10 times
    unsigned int max_iterations = 10000000; // repeat at most 1e7 times
    Timer timer; // timer to collect runtime of each individual run

    // vector that stores the timing of each benchmark case
    std::vector<accumulator_t> timings_eigen;
    std::vector<accumulator_t> timings_strassen;

    // store dimension of the matrix
    std::vector<unsigned> matrix_dims;

    // display header column
    std::cout << std::setw(4) << "2^k" 
              << std::setw(15) << "A*B" 
              << std::setw(15) << "Strassen" << std::endl;

    for(unsigned k = 4; k <= 9; k++) {
        // initialize new accumulator
        accumulator_t time_accumulator_eigen,
                      time_accumulator_strassen;

        unsigned int n = pow(2,k);
        matrix_dims.push_back(n);

        // benchmark with a matrix size of 2^k for `min_runtime` seconds
        {
            unsigned it = 0; // iteration counter
            srand(seed); // seed random number generator
            // repeat benchmark for `min_runtime` with the constraint that
            //  the number of iterations is less then max_iterations but 
            //  at least `min_iterations`
            while (it < max_iterations
                    && (it < min_iterations
                        || sum(time_accumulator_eigen) < min_runtime)) {
                // initialize random input matricies
                MatrixXd A = MatrixXd::Random(n, n);
                MatrixXd B = MatrixXd::Random(n, n);

                // initialize memory for result matrix
                MatrixXd AxB(n,n);
                
                // benchmark eigens matrix multiplication
                timer.start(); // start timer
                AxB=A*B; // do the multiplication
                escape(AxB); // make sure that the calculation is not optimized away
                timer.stop(); // stop timer
                time_accumulator_eigen(timer.duration()); // store timing

                ++it; // step iteration counter
            }
            
            // store accumulator for later plotting
            timings_eigen.push_back(time_accumulator_eigen);
        }
        {
            unsigned it = 0; // iteration counter
            srand(seed); // seed random number generator
            // repeat benchmark for `min_runtime` with the constraint that
            //  the number of iterations is less then max_iterations but 
            //  at least `min_iterations`
            while (it < max_iterations
                    && (it < min_iterations
                        || sum(time_accumulator_strassen) < min_runtime)) {
                // initialize random input matricies
                MatrixXd A = MatrixXd::Random(n, n);
                MatrixXd B = MatrixXd::Random(n, n);

                // initialize memory for result matrix
                MatrixXd AxB(n,n);

                // benchmark stassens matrix multiplication
                timer.start(); // start timer
                AxB=strassenMatMult(A, B); // do the multiplication
                escape(AxB); // make sure that the calculation is not optimized away
                timer.stop(); // stop timer
                time_accumulator_strassen(timer.duration()); // store timing

                ++it; // step iteration counter
            }

            // store accumulator for later plotting
            timings_strassen.push_back(time_accumulator_strassen);
        }

        // print mean runtime
        std::cout << std::setw(4) << k // power
                  << std::setw(15) << mean(time_accumulator_eigen) << std::scientific // eigen timing
                  << std::setw(15) << mean(time_accumulator_strassen) << std::scientific // strassing timing
                  << std::endl;
    }

    // calculate asymptotic complexity
    std::cout << "Eigen: O(n^" << asymptotic_complexity(matrix_dims, timings_eigen, 2)[0] << ")" << std::endl;
    std::cout << "Strassen: O(n^" << asymptotic_complexity(matrix_dims, timings_strassen, 2)[0] << ")" << std::endl;
}
