//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Build an "arrow matrix" and compute A*A*y
 * Given vectors $a$ and $d$, returns A*A*x in $y$, where A is built from a, d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
void arrow_matrix_2_times_x(const VectorXd &d, const VectorXd &a,
                            const VectorXd &x, VectorXd &y) {
    assert(d.size() == a.size() && a.size() == x.size() &&
           "Vector size must be the same!");
    int n = d.size();

    VectorXd d_head = d.head(n-1);
    VectorXd a_head = a.head(n-1);
    MatrixXd d_diag = d_head.asDiagonal();

    MatrixXd A(n,n);

    A << d_diag,             a_head,
         a_head.transpose(), d(n-1);

    y = A*A*x;
}

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d An n-dimensional vector
 * @param[in] a An n-dimensional vector
 * @param[in] x An n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
void efficient_arrow_matrix_2_times_x(const VectorXd &d,
                                      const VectorXd &a,
                                      const VectorXd &x,
                                      VectorXd &y) {
    assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
    int n = d.size();

    // TODO: Implement an efficient version of arrow\_matrix\_2\_times\_x
}

/* \brief Compute the runtime of arrow matrix multiplication.
 * Repeat tests 10 times, and output the minimal runtime
 * amongst all times. Test both the inefficient and the efficient
 * versions.
*/
void runtime_arrow_matrix() {
    // TODO: your code here, time the codes

}


int main(void) {
    // Test vectors
    VectorXd a(5);
    a << 1., 2., 3., 4., 5.;
    VectorXd d(5);
    d <<1., 3., 4., 5., 6.;
    VectorXd x(5);
    x << -5., 4., 6., -8., 5.;
    VectorXd yi;

    // Run both functions
    arrow_matrix_2_times_x(a,d,x,yi);
    VectorXd ye(yi.size());
    efficient_arrow_matrix_2_times_x(a,d,x,ye);

    // Compute error
    double err = (yi - ye).norm();

    // Output error
    std::cout << "--> Correctness test." << std::endl;
    std::cout << "Error: " << err << std::endl;

    // Print out runtime
    std::cout << "--> Runtime test." << std::endl;
    runtime_arrow_matrix();

    // Final test: exit with error if error is too big
    double eps = std::numeric_limits<double>::denorm_min();
    exit(err < eps);
}
