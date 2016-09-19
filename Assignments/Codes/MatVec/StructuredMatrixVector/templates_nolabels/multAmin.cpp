//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <Eigen/Dense>

#include <iostream>
#include <vector>

#include "timer.h"

using namespace Eigen;

//! \brief build A*x using array of ranges and of ones.
//! \param[in] x vector x for A*x = y
//! \param[out] y y = A*x
void multAminSlow(const Eigen::VectorXd & x, Eigen::VectorXd & y) {
    unsigned int n = x.size();

    VectorXd one = Eigen::VectorXd::Ones(n);
    VectorXd linsp = VectorXd::LinSpaced(n,1,n);
    y = ( (one * linsp.transpose()).cwiseMin(linsp * one.transpose()) ) *x;
}

//! \brief build A*x using a simple for loop
//! \param[in] x vector x for A*x = y
//! \param[out] y y = A*x
void multAminLoops(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();

    MatrixXd A(n,n);

    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            A(i,j) = std::min(i+1,j+1);
        }
    }
    y = A * x;
}

//! \brief build A*x using a clever representation
//! \param[in] x vector x for A*x = y
//! \param[out] y y = A*x
void multAmin(const VectorXd & x, VectorXd & y) {
    unsigned int n = x.size();
    y = VectorXd::Zero(n);
    VectorXd v = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);

    v(0) = x(n-1);
    w(0) = x(0);

    for(unsigned int j = 1; j < n; ++j) {
        v(j) = v(j-1) + x(n-j-1);
        w(j) = w(j-1) + (j+1)*x(j);
    }
    for(unsigned int j = 0; j < n-1; ++j) {
        y(j) = w(j) + v(n-j-2)*(j+1);
    }
    y(n-1) = w(n-1);
}

int main(void) {
    // Build Matrix B with 10x10 dimensions such that B = inv(A)
    unsigned int n = 10;
    MatrixXd B = MatrixXd::Zero(n,n);
    for(unsigned int i = 0; i < n; ++i) {
        B(i,i) = 2;
        if(i < n-1) B(i+1,i) = -1;
        if(i > 0) B(i-1,i) = -1;
    }
    B(n-1,n-1) = 1;
    std::cout << "B = " << B << std::endl;

    // Check that B = inv(A) (up to machine precision)
    VectorXd x = VectorXd::Random(n), y;
    multAmin(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminSlow(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;
    multAminLoops(B*x, y);
    std::cout << "|y-x| = " << (y - x).norm() << std::endl;

    // Timing from 2^4 to 2^13 repeating nruns times
    unsigned int nruns = 10;
    std::vector<double> times_slow, times_slow_loops, times_fast;
    for(unsigned int p = 4; p <= 13; ++p) {
        Timer tm_slow, tm_slow_loops, tm_fast;
        for(unsigned int r = 0; r < nruns; ++r) {
            x = VectorXd::Random(pow(2,p));

            tm_slow.start();
            multAminSlow(x, y);
            tm_slow.stop();

            tm_slow_loops.start();
            multAminLoops(x, y);
            tm_slow_loops.stop();

            tm_fast.start();
            multAmin(x, y);
            tm_fast.stop();
        }
        times_slow.push_back( tm_slow.mean() );
        times_slow_loops.push_back( tm_slow_loops.mean() );
        times_fast.push_back( tm_fast.mean() );
    }

    for(auto it = times_slow.begin(); it != times_slow.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_slow_loops.begin(); it != times_slow_loops.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    for(auto it = times_fast.begin(); it != times_fast.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
}
