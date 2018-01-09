//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "timer.h"

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Efficiently compute the Matrix-Vector-Product
 * $x=Ay$, where A(i,j) is zero except for values on the
 * diagonal and antidiagonal. There the values are given
 * by the entries of the vector a.
 *
 * @param[in]  a, vector storing the elements of $A$
 * @param[in]  y, vector to multiply $A$ with
 * @param[out] x, result vector of $Ay=x$
 */
template<class Vector>
void xmatmult(const Vector& a, const Vector& y, Vector& x) {
    assert(a.size() == y.size() && a.size() == x.size()
            && "Input vector dimensions must match");
    unsigned n = a.size();
    
    // TODO: Efficiently compute A*y
    // and store the result in x
}

/*@brief Compares the runtimes of the efficient
 * multipication to the normal Matrix-Vector product
 */
void compare_times() {
    std::cout << "Measuring runtimes for comparison" << std::endl;
    // TODO: Measure the runtime of your implementation and that of Eigen's
    // Matrix-Vector multiplication with the timer.h library for
    // n = 2^5,2^6,...2^14

}

void test() {
    // testing for even n
    unsigned n = 10;
    VectorXd a,y,x(n);
    a = y = VectorXd::Random(n,1);
    //building A for normal Matrix-Vector multiplication O(n*n)
    MatrixXd A = a.asDiagonal();
    for (unsigned i = 0; i < n; ++i) {
        A(n-i-1,i) = A(i,i);
    }
    xmatmult(a,y,x);
    
    double error = (x - A*y).norm();
    if (std::abs(error) > 1e-10) {    //== on floating point numbers is meaninless
        std::cout << "Wrong result for even n" << std::endl;
        return;
    }
    
    // testing for odd n
    n = 11;
    a = y = x= VectorXd::Random(n,1);
    A = a.asDiagonal();
    for (unsigned i = 0; i < n; ++i) {
        A(n-i-1,i) = A(i,i);
    }
    
    xmatmult(a,y,x);
    
    error = (x - A*y).norm();
    if (std::abs(error) > 1e-10) {
        std::cout << "Wrong result for odd n" << std::endl;
        return;
    }
    
    std::cout << "No errors found in your implementation" << std::endl;
}

int main() {
    test();
    compare_times();
}
