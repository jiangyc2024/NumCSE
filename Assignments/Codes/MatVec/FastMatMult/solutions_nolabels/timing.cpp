//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

#include "strassen.hpp"

int main() {
    // Store seed for rng
    unsigned seed = (unsigned) time(0);
     // seed random number generator
    srand(seed);

    // Minimum number of repetitions
    unsigned int repetitions = 10;


    // Display header column
    std::cout << std::setw(4)  << "k"
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

        // Print runtimes
        std::cout << std::setw(4) << k // power
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << timer.min() // eigen timing
                  << std::setprecision(3) << std::setw(15) << std::scientific
                  << timer_own.min() // strassing timing
                  << std::endl;

    }

}
