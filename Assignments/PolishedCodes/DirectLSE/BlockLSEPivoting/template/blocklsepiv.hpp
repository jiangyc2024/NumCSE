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

  // TODO: (2-2.b) Efficiently calculate A*x. Eigen's array arithmetic is
  // helpful here.
  // START

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
  Eigen::VectorXd x = Eigen::VectorXd::Zero(2 * n);
  constexpr double eps = std::numeric_limits<double>::epsilon();

  // TODO: (2-2.d) Solve Ax = b for x using Gaussian elimination without
  // directly envoking Eigen's solvers on the full matrix
  // START

  // END

  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Tabulates runtimes for different n and plots them.
 *
 */
/* SAM_LISTING_BEGIN_3 */
void numericalExperiment() {
  constexpr unsigned int repeats = 5;
  Eigen::VectorXd d1, d2, c, b, y;
  plt::figure();

  // TODO: (2-2.f) Tabulate runtimes of solveA and plot them using
  // matplotlibcpp.
  // START

  // END
  plt::savefig("cx_out/complexity.png");
}
/* SAM_LISTING_END_3 */

#endif
