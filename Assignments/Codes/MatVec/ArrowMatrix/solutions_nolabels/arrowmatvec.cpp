//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
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

    // In this lines, we extract the blocks used to construct the matrix A.
    VectorXd d_head = d.head(n-1);
    VectorXd a_head = a.head(n-1);
    MatrixXd d_diag = d_head.asDiagonal();

    MatrixXd A(n,n);

    // We build the matrix A using the "comma initialization": each expression separated
    // by a comma is a "block" of the matrix we are building.
    // d\_diag is the top left (n-1)x(n-1) block
    // a\_head is the top right vertical vector
    // a\_head.transpose() is the bottom left horizontal vector
    // d(n-1) is a single element (a 1x1 matrix), on the bottom right corner
    // This is how the matrix looks like:
    // A = | D   | a      |
    //     |-----+--------|
    //     | a\^T | d(n-1) |
    A << d_diag,             a_head,
         a_head.transpose(), d(n-1);

    y = A*A*x;
}

/* @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 * @param[in] d A n-dimensional vector
 * @param[in] a A n-dimensional vector
 * @param[in] x A n-dimensional vector
 * @param[out] y The vector y = A*A*x
 */
void efficient_arrow_matrix_2_times_x(const VectorXd &d,
                                      const VectorXd &a,
                                      const VectorXd &x,
                                      VectorXd &y) {
    assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
    int n = d.size();

    // Notice that we can compute (A*A)*x more efficiently using
    // A*(A*x). This is, in fact, performing two matrix vector
    // multiplications
    // instead of a more expensive matrix-matrix multiplication.
    // Therefore, as first step, we need a way to efficiently
    // compute A*x

    // This function computes A*x. you can use it
    // by calling A\_times\_x(x).
    // This is the syntax for lambda functions: notice the extra
    // [variables] code. Each variable written within [] brackets
    // will be captured (i.e. seen) inside the lambda function.
    // Without \&, a copy of the variable will be performed.
    // Notice that A = D + H + V, s.t. A*x = D*x + H*x + V*x
    // D*x can be rewritten as d*x componentwise
    // H*x is zero, except at the last component
    // V*x is only affected by the last component of x
    auto A_times_x = [&a, &d, n] (const VectorXd & x) {
        // This takes care of the diagonal (D*x)
        // Notice: we use d.array() to tell Eigen to treat
        // a vector as an array. As a result: each operation
        // is performed componentwise.
        VectorXd Ax = ( d.array() * x.array() ).matrix();

        // H*x only affects the last component of A*x
        // This is a dot product between a and x with the last
        // component removed
        Ax(n-1) += a.head(n - 1).dot(x.head(n-1));

        // V*x is equal to the vector
        // (a(0)*x(n-1), ..., a(n-2)*x(n-1), 0)
        Ax.head(n-1) +=  x(n-1) * a.head(n-1);

        return Ax;
    };

    // <=> y = A*(A*x)
    y = A_times_x(A_times_x(x));
}

/* \brief Compute the runtime of arrow matrix multiplication.
 * Repeat tests 10 times, and output the minimal runtime
 * amongst all times. Test both the inefficient and the efficient
 * versions.
*/
void runtime_arrow_matrix() {

    std::cout << std::setw(8) << "n"
              << std::setw(15) << "original"
              << std::setw(15) << "efficient"
              << std::endl;
    for (unsigned n = 4; n <= 2048; n *= 2) {
        // Number of repetitions
        unsigned int repeats = 10;

        Timer timer, timer_eff;
        // Repeat test many times
        for (unsigned int r = 0; r < repeats; ++r) {
            // Create test input using random vectors
            Eigen::VectorXd a = Eigen::VectorXd::Random(n),
                            d = Eigen::VectorXd::Random(n),
                            x = Eigen::VectorXd::Random(n),
                            y;

            // Compute times
            timer.start();
            arrow_matrix_2_times_x(d, a, x, y);
            timer.stop();

            // Compute times for efficient implementation
            timer_eff.start();
            efficient_arrow_matrix_2_times_x(d, a, x, y);
            timer_eff.stop();
        }


        // Print results
        std::cout << std::setw(8) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << timer.min()
                  << std::setw(15) << timer_eff.min()
                  << std::endl;
    }

}

// Rename long variable name to duration_t (easy to change)
using duration_t = std::chrono::nanoseconds;

/* \brief Compute runtime of $F$, repeating test "repeats" times
 * Will return minimal runtime.
 * This function uses "crhono".
 * \tparam Function type of F, must have an operator()
 * \param[in] F Function for which you want to measure runtime.
 * \param[in] repeats Number of repetitions.
 */
template <class Function>
duration_t timing(const Function & F, int repeats = 10) {
    // Shortcut for time_point
    using time_point_t = std::chrono::high_resolution_clock::time_point;

    // Loop many times
    duration_t min_elapsed;
    for(int r = 0; r < repeats; r++) {
        // Start clock (MATLAB: tic)
        time_point_t start = std::chrono::high_resolution_clock::now();

        // Run function
        F();

        // Stop clock (MATLAB: toc) and measure difference
        duration_t elapsed = std::chrono::duration_cast<duration_t>(
                            std::chrono::high_resolution_clock::now() - start);

        // Compute min between all runs
        min_elapsed = r == 0 ? elapsed : std::min(elapsed, min_elapsed);
    }

    return min_elapsed;

}

/* \brief Compute timing using chrono
 * Also demonstrate use of lambda functions
 */
void runtime_arrow_matrix_with_chrono() {

    // Table header
    std::cout << std::setw(8) << "n"
              << std::scientific << std::setprecision(3)
              << std::setw(15) << "original"
              << std::setw(15) << "efficient"
              << std::endl;

    // Run from $2^5$ to $2^11$ with powers of two
    for(unsigned int n = (1 << 5); n < (1 << 12); n = n << 1) {
        // Create random vectors
        VectorXd d = VectorXd::Random(n);
        VectorXd a = VectorXd::Random(n);
        VectorXd x = VectorXd::Random(n);
        VectorXd y(n);

        // Call "timing", using a lambda function for F
        // Remember: we cannot pass arrow\_matrix\_2\_times\_x directly to timing
        // the timing function expects a n object with operator()(void)
        duration_t elapsed = timing([&a, &d, &x, &y] () {
            arrow_matrix_2_times_x(d, a, x, y);
        }, 10);
        // Call "timing", using a lambda function for F
        duration_t elapsed_efficient = timing([&a, &d, &x, &y] () {
            efficient_arrow_matrix_2_times_x(d, a, x, y);
        }, 10);

        // Output timings
        std::cout << std::setw(8)<< n
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << elapsed.count() * 1e-9 // ns to s
                  << std::setw(15) << elapsed_efficient.count() * 1e-9 // ns to s
                  << std::endl;
    }

}

int main(void) {
    // Test vectors
    VectorXd a(5);
    a << 1., 2., 3., 4., 5.;
    VectorXd d(5);
    d <<1., 3., 4., 5., 6.;
    VectorXd x(5);
    x << -5., 4., 6., -8., 5.;
    VectorXd yi, ye;

    // Run both functions
    arrow_matrix_2_times_x(a,d,x,yi);
    efficient_arrow_matrix_2_times_x(a,d,x,ye);

    // Compute error
    double err = (yi - ye).norm();

    // Output error
    std::cout << "--> Correctness test." << std::endl;
    std::cout << "Error: " << err << std::endl;

    // Print out runtime
    std::cout << "--> Runtime test." << std::endl;
    runtime_arrow_matrix();
    // Print out runtime with chrono
    std::cout << "--> Runtime test. with chrono." << std::endl;
    runtime_arrow_matrix_with_chrono();

    // Final test: exit with error if error is too big
    double eps = std::numeric_limits<double>::denorm_min();
    exit(err < eps);
}
