#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

// boost accumulators
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "timer.h"

#include "strassen.hpp"
#include "polyfit.hpp"

using namespace Eigen;

// we use accumulators to collect important characteristics
//  like minimum-, maximum- and mean-runtime
using namespace boost::accumulators;

using accumulator_t = boost::accumulators::accumulator_set<
    double, // collect number of seconds in double precision
    stats<tag::mean> // collect mean
>;

// tell the compiler that `p` is used in some sense and should not be
//  optimized away
template <typename T>
static void escape (T&& p) {
    asm volatile("" :  : "g"(p) : "memory");
}

// estimate asymptotic computational complexity
VectorXd asymptotic_complexity(std::vector<unsigned> dims, std::vector<double>& means,
  unsigned offset=0) {
    assert(dims.size() == means.size());
    assert(dims.size() > offset);
    // convert std::vector<unsigned> into Eigen::VectorXd
    VectorXd means_ = Map<Matrix<double, 1, -1>>(means.data()+offset, 1,
        means.size()-offset);
    VectorXd dims_ = Map<Matrix<unsigned, 1, -1>>(dims.data()+offset, 1,
        dims.size()-offset).cast<double>();

    return polyfit(dims_.array().log().matrix(), means_.array().log().matrix(), 1);
}

int main()
{
    unsigned seed = (unsigned) time(0); // store seed for rng
    srand(seed); // seed random number generator

    double min_runtime = 20.; // minimum runtime in seconds
    unsigned int min_iterations = 10; // repeat at least 10 times
    unsigned int max_iterations = 10000000; // repeat at most 1e7 times
    Timer timer; // timer to collect runtime of each individual run

    // vector that stores the timing of each benchmark case
    std::vector<double> timings_eigen;
    std::vector<double> timings_strassen;

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

        // initialize random input matricies
        MatrixXd A = MatrixXd::Random(n, n);
        MatrixXd B = MatrixXd::Random(n, n);

        // benchmark with a matrix size of 2^k for `min_runtime` seconds
        {
            unsigned it = 0; // iteration counter

            // repeat benchmark for `min_runtime` with the constraint that
            //  the number of iterations is less then max_iterations but
            //  at least `min_iterations`
            while (it < max_iterations
                    && (it < min_iterations
                        || sum(time_accumulator_eigen) < min_runtime)) {
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
            timings_eigen.push_back(mean(time_accumulator_eigen));
        }
        {
            unsigned it = 0; // iteration counter

            // repeat benchmark for `min_runtime` with the constraint that
            //  the number of iterations is less then max_iterations but
            //  at least `min_iterations`
            while (it < max_iterations
                    && (it < min_iterations
                        || sum(time_accumulator_strassen) < min_runtime)) {
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
            timings_strassen.push_back(mean(time_accumulator_strassen));
        }

        // print mean runtime
        std::cout << std::setw(4) << k // power
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << mean(time_accumulator_eigen) // eigen timing
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << mean(time_accumulator_strassen) // strassing timing
                  << std::endl;
    }

    // estimate asymptotic computational complexity
    std::cout << "Estimated asymptotic complexity:" << std::endl;
    std::cout << "Eigen: O(n^" << asymptotic_complexity(matrix_dims, timings_eigen, 2)[0] << ")" << std::endl;
    std::cout << "Strassen: O(n^" << asymptotic_complexity(matrix_dims, timings_strassen, 2)[0] << ")" << std::endl;
}
