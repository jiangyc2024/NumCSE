#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

#include "strassen.hpp"

int main() {
    /* SAM_LISTING_BEGIN_1 */
    // Store seed for rng
    unsigned seed = (unsigned) time(0);
     // seed random number generator
    srand(seed);

    // Minimum number of repetitions
    unsigned int repetitions = 10;

#if SOLUTION
    // Display header column
    std::cout << std::setw(4)  << "2^k"
              << std::setw(15) << "A*B"
              << std::setw(15) << "Strassen" << std::endl;
    for(unsigned k = 4; k <= 9; k++) {
        unsigned int n = std::pow(2, k);

        // Initialize random input matricies
        MatrixXd A = MatrixXd::Random(n, n);
        MatrixXd B = MatrixXd::Random(n, n);

        // Timer to collect runtime of each individual run
        Timer timer, timer_own;

        for(unsigned int r = 0; r < repetitions; ++r) {
            // initialize memory for result matrix
            MatrixXd AxB, AxB2;

            // Benchmark eigens matrix multiplication
            timer.start(); // start timer
            AxB=(A*B); // do the multiplication
            timer.stop(); // stop timer

            // Benchmark Stassens matrix multiplication
            timer_own.start(); // start timer
            AxB2=strassenMatMult(A, B); // do the multiplication
            timer_own.stop(); // stop timer

            // volatile double a = (AxB+AxB2).norm();
        }

        // Print runtime
        std::cout << std::setw(4) << k // power
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << timer.min() // eigen timing
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << timer_own.min() // strassing timing
                  << std::endl;
    }
#else // TEMPLATE
    // TODO: time Strassen and Eigen matrix multiplication
    // repeat 10 times and output the minimum runtime
    // use 3 digits and scientific notation
#endif // TEMPLATE
    /* SAM_LISTING_END_1 */
}
