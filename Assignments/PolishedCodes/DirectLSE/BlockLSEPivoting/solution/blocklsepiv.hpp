#ifndef BLOCKLSEPIV_HPP
#define BLOCKLSEPIV_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "matplotlibcpp.h"
#include "timer.h"

namespace plt = matplotlibcpp;

/**
 * @brief Multiplies A*x. Build matrix with block structure.
 *
 * @param d1 n-dim vector
 * @param d2 n-dim vector
 * @param c n-dim vector
 * @param x vector with dim 2n
 * @return Eigen::VectorXd A*x
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd multA(const Eigen::VectorXd& d1, const Eigen::VectorXd& d2,
                      const Eigen::VectorXd& c, const Eigen::VectorXd& x) {
  const unsigned int n = d1.size();
  assert(n == d2.size() && n == c.size() && 2 * n == x.size() &&
         "Size mismatch!");
  Eigen::ArrayXd y = Eigen::VectorXd::Zero(2 * n);

  // TODO: (3-2.b) Efficiently calculate A*x. Eigen's array arithmetic is
  // helpful here.
  // START
  auto& x1 = x.head(n).array();
  auto& x2 = x.tail(n).array();

  y << d1.array() * x1 + c.array() * x2, d2.array() * x2 + c.array() * x1;
  // END
  return y;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves A*x = b for x using partial pivoting and exploiting the matrix'
 * structure.
 *
 * @param d1 n-dim vector
 * @param d2 n-dim vector
 * @param c n-dim vector
 * @param b r.h.s. vector with dim 2n
 * @return Eigen::VectorXd x s.t. A*x = b
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveA(const Eigen::VectorXd& d1, const Eigen::VectorXd& d2,
                       const Eigen::VectorXd& c, const Eigen::VectorXd& b) {
  const unsigned int n = d1.size();
  assert(n == d2.size() && n == c.size() && 2 * n == b.size() &&
         "Size mismatch!");

  Eigen::ArrayXd c1_arr = c, c2_arr = c, d1_arr = d1, d2_arr = d2, b_arr = b;

  constexpr double eps = std::numeric_limits<double>::epsilon();

  // TODO: (3-2.d) Solve Ax = b for x using Gaussian elimination without using
  // Eigen's built-in solvers
  // START

  // for forward elimination and pivoting we are only required to loop the first
  // half of the matrix; loop over diagonal
  for (unsigned int k = 0; k < n; ++k) {
    // check if need to pivot (i.e. swap two rows)
    // luckily we only need to check two rows
    const double maxk = std::max(std::abs(d1_arr(k)), std::abs(c1_arr(k)));
    const double maxnk = std::max(std::abs(c2_arr(k)), std::abs(d2_arr(k)));
    if (std::abs(c1_arr(k)) / maxk > std::abs(d1_arr(k)) / maxnk) {
      // matrix
      std::swap(d1_arr(k), c2_arr(k));
      std::swap(c1_arr(k), d2_arr(k));
      // r.h.s.
      std::swap(b_arr(n), b_arr(n + k));
    }

    // check if matrix is almost singular
    const double piv = d1_arr(k);
    // norm of the block from k, k to n - 1, n - 1
    const double norm = std::abs(d1_arr(k)) + std::abs(d2_arr(k)) +
                        std::abs(c1_arr(k)) + std::abs(c2_arr(k));
    if (piv < eps * norm) {
      std::cout << "Warning: matrix nearly singular!" << std::endl;
    }

    // Multiplication factor
    const double fac = c2_arr(k) / piv;

    // Actually perform substitution
    // botton right portion changes
    d2_arr(k) -= c1_arr(k) * fac;
    // r.h.s.
    b_arr(n + k) -= b_arr(k) * fac;
  }

  // Now the system has the form:
  // | d1 | c  |   |   |   |   |
  // | 0  | d2 | * | x | = | b |
  // with d1, d2, c diagonal

  // Backward substitution
  // lower portion
  b_arr.tail(n) /= d2_arr;
  // upper portion
  b_arr.head(n) = (b_arr.head(n) - c1_arr * b_arr.tail(n)) / d1_arr;
  // END

  return b_arr;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Tabulates runtimes for different n and plots them.
 *
 */
/* SAM_LISTING_BEGIN_3 */
void numericalExperiment() {
  constexpr unsigned int repeats = 10;
  Eigen::VectorXd d1, d2, c, b, y;
  plt::figure();

  // TODO: (3-2.f) Tabulate runtimes of solveA and plot them using
  // matplotlibcpp. START
  std::vector<double> n(7), times(7), o_n(7);

  // Table header
  std::cout << std::setw(10) << "k" << std::scientific << std::setprecision(3)
            << std::setw(15) << "runtime [s]" << std::endl;

  // Loop from $2$ to $2^7$
  for (unsigned int k = 2, i = 0; k <= (1 << 7); k <<= 1, ++i) {
    d1 = Eigen::VectorXd::LinSpaced(k, 1, k);
    d2 = -d1;
    c = Eigen::VectorXd::Ones(k);
    b.resize(2 * k);
    b << d1, d1;

    // Repeat measurements
    Timer tm;
    for (unsigned int r = 0; r < repeats; ++r) {
      // Test run
      tm.start();
      y = solveA(d1, d2, c, b);
      tm.stop();
    }

    // Keep track of data for plotting
    n[i] = k;
    times[i] = tm.duration();
    o_n[i] = k * 1e-5;  // scale to beautify plot

    // Print measurements
    std::cout << std::setw(10) << k << std::scientific << std::setprecision(3)
              << std::setw(15) << tm.duration() << std::endl;
  }

  // Plotting
  plt::loglog(n, times, "o", {{"label", "solveA"}});
  plt::loglog(n, o_n, {{"label", "O(n)"}});
  plt::title("timings for solveA accumulated over 10 runs");
  plt::xlabel("dimension n");
  plt::ylabel("time [s]");
  plt::legend();
  // END
  plt::savefig("cx_out/complexity.png");
}
/* SAM_LISTING_END_3 */

#endif