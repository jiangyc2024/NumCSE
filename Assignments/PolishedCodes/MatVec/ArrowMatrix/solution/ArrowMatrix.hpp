#ifndef ARROWMATRIX_HPP
#define ARROWMATRIX_HPP

#include <Eigen/Dense>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "plot.hpp"
#include "timer.h"

/**
 * @brief Build an "arrow matrix" and compute A*A*y
 * Given vectors $a$ and $d$, returns A*A*x in $y$, where A is built from a, d
 *
 * @param d An n-dimensional vector
 * @param a An n-dimensional vector
 * @param x An n-dimensional vector
 * @param y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_0 */
void arrow_matrix_2_times_x(const Eigen::VectorXd &d, const Eigen::VectorXd &a,
                            const Eigen::VectorXd &x, Eigen::VectorXd &y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  const unsigned int n = d.size();

  // In this lines, we extract the blocks used to construct the matrix A.
  Eigen::VectorXd d_head = d.head(n - 1);
  Eigen::VectorXd a_head = a.head(n - 1);
  Eigen::MatrixXd d_diag = d_head.asDiagonal();

  Eigen::MatrixXd A(n, n);

  // We build the matrix A using the "comma initialization": each expression
  // separated by a comma is a "block" of the matrix we are building. d\_diag is
  // the top left (n-1)x(n-1) block, a\_head is the top right vertical vector,
  // a\_head.transpose() is the bottom left horizontal vector,
  // d(n-1) is a single element (a 1x1 matrix) on the bottom right corner.
  // This is how the matrix looks like:
  // A = | D    | a      |
  //     |------+--------|
  //     | a\^T | d(n-1) |
  A << d_diag, a_head, a_head.transpose(), d(n - 1);

  y = A * A * x;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Build an "arrow matrix"
 * Given vectors $a$ and $b$, returns A*A*x in $y$, where A is build from a,d
 *
 * @param d An n-dimensional vector
 * @param a An n-dimensional vector
 * @param x An n-dimensional vector
 * @param y The vector y = A*A*x
 */
/* SAM_LISTING_BEGIN_1 */
void efficient_arrow_matrix_2_times_x(const Eigen::VectorXd &d,
                                      const Eigen::VectorXd &a,
                                      const Eigen::VectorXd &x,
                                      Eigen::VectorXd &y) {
  assert(d.size() == a.size() && a.size() == x.size() &&
         "Vector size must be the same!");
  const unsigned int n = d.size();

  // TODO: (1-1.c) Implement an efficient version of
  // arrow\_matrix\_2\_times\_x. Hint: Notice that performing two matrix vector
  // multiplications is less expensive than a matrix-matrix multiplication.
  // Therefore, as first step, we need a way to efficiently compute A*x

  // START

  // Notice that we can compute (A*A)*x more efficiently using
  // A*(A*x). This is, in fact, performing two matrix vector
  // multiplications instead of a more expensive matrix-matrix multiplication.
  // Therefore, as a first step, we need a way to efficiently
  // compute A*x

  // This function computes A*x. You can use it by calling A\_times\_x(x).
  // This is the syntax for lambda functions: notice the extra
  // [variables] code. Each variable written within [] brackets
  // will be captured (i.e. seen) inside the lambda function.
  // Without \&, a copy of the variable will be performed.
  // Notice that A = D + H + V, s.t. A*x = D*x + H*x + V*x
  // D*x can be rewritten as d*x componentwise
  // H*x is zero, except at the last component
  // V*x is only affected by the last component of x
  auto A_times_x = [&a, &d, n](const Eigen::VectorXd &x) {
    // This takes care of the diagonal (D*x)
    // Notice: we use d.array() to tell Eigen to treat
    // a vector as an array. As a result: each operation
    // is performed componentwise.
    Eigen::VectorXd Ax = (d.array() * x.array()).matrix();

    // H*x only affects the last component of A*x
    // This is a dot product between a and x with the last
    // component removed
    Ax(n - 1) += a.head(n - 1).dot(x.head(n - 1));

    // V*x is equal to the vector
    // (a(0)*x(n-1), ..., a(n-2)*x(n-1), 0)
    Ax.head(n - 1) += x(n - 1) * a.head(n - 1);

    return Ax;
  };

  // <=> y = A*(A*x)
  y = A_times_x(A_times_x(x));

  // END
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute the runtime of arrow matrix multiplication.
 * Repeat tests 10 times, and output the minimal runtime
 * amongst all times. Test both the inefficient and the efficient
 * versions.
 *
 */
/* SAM_LISTING_BEGIN_3 */
void runtime_arrow_matrix() {
  // Memory allocation for plot
  std::vector<double> vec_size;
  std::vector<double> elap_time, elap_time_eff;

  // header for result print out
  std::cout << std::setw(8) << "n" << std::setw(15) << "original"
            << std::setw(15) << "efficient" << std::endl;

  for (unsigned int n = 32; n <= 1024; n <<= 1) {
    // save vector size (for plot)
    vec_size.push_back(n);

    // Number of repetitions
    constexpr unsigned int repeats = 10;

    Timer timer, timer_eff;
    // Repeat test many times
    for (unsigned int r = 0; r < repeats; ++r) {
      // Create test input using random vectors
      Eigen::VectorXd a = Eigen::VectorXd::Random(n);
      Eigen::VectorXd d = Eigen::VectorXd::Random(n);
      Eigen::VectorXd x = Eigen::VectorXd::Random(n);
      Eigen::VectorXd y;

      // Compute times for original implementation
      timer.start();
      arrow_matrix_2_times_x(d, a, x, y);
      timer.stop();

      // TODO: (1-1.e) Compute times for efficient implementation
      // START
      timer_eff.start();
      efficient_arrow_matrix_2_times_x(d, a, x, y);
      timer_eff.stop();
      // END
    }

    // Print results (for grading): inefficient
    std::cout << std::setw(8) << n << std::scientific << std::setprecision(3)
              << std::setw(15) << timer.min() << std::setw(15)
              << timer_eff.min() << std::endl;

    // time needed
    elap_time.push_back(timer.min());
    elap_time_eff.push_back(timer_eff.min());
  }

  /* DO NOT CHANGE */
  // create plot
  plot(vec_size, elap_time, elap_time_eff, "./cx_out/comparison.png");
}
/* SAM_LISTING_END_3 */

// Additional

// Rename long variable name to duration_t (easy to change)
using duration_t = std::chrono::nanoseconds;

/**
 * @brief Compute runtime of $F$, repeating test "repeats" times
 * Will return minimal runtime. This function uses "chrono".
 *
 * @tparam Function type of F, must have an operator()
 * @param F Function for which you want to measure runtime.
 * @param repeats Number of repetitions.
 * @return duration_t runtime of F
 */
template <class Function>
duration_t timing(const Function &F, unsigned int repeats = 10) {
  // Shortcut for time_point
  using time_point_t = std::chrono::high_resolution_clock::time_point;

  // Loop many times
  duration_t min_elapsed;
  for (unsigned int r = 0; r < repeats; r++) {
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

/**
 * @brief Compute timing using chrono
 * Also demonstrate use of lambda functions
 *
 */
void runtime_arrow_matrix_with_chrono() {
  // Table header
  std::cout << std::setw(8) << "n" << std::scientific << std::setprecision(3)
            << std::setw(15) << "original" << std::setw(15) << "efficient"
            << std::endl;

  // Run from $2^3$ to $2^7$ with powers of two
  for (unsigned int n = 8; n <= 128; n <<= 1) {
    // Create random vectors
    Eigen::VectorXd d = Eigen::VectorXd::Random(n);
    Eigen::VectorXd a = Eigen::VectorXd::Random(n);
    Eigen::VectorXd x = Eigen::VectorXd::Random(n);
    Eigen::VectorXd y(n);

    // Call "timing", using a lambda function for F
    // Remember: we cannot pass arrow\_matrix\_2\_times\_x directly to timing
    // the timing function expects a n object with operator()(void)
    duration_t elapsed =
        timing([&a, &d, &x, &y]() { arrow_matrix_2_times_x(d, a, x, y); }, 10);
    // Call "timing", using a lambda function for F
    duration_t elapsed_efficient = timing(
        [&a, &d, &x, &y]() { efficient_arrow_matrix_2_times_x(d, a, x, y); },
        10);

    // Output timings (not part of the exercice)
    std::cout << std::setw(8) << n << std::scientific << std::setprecision(3)
              << std::setw(15) << elapsed.count() * 1e-9            // ns to s
              << std::setw(15) << elapsed_efficient.count() * 1e-9  // ns to s
              << std::endl;
  }
}

#endif
