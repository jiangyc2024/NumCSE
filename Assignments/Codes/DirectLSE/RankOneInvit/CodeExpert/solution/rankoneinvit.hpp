#ifndef RANKONEINVIT_HPP
#define RANKONEINVIT_HPP

// TO DO: (3-10.a) Write the function rankoneinvit().
// (3-10.c) Write the function rankoneinvit_fast().
// (3-10.e) Write a function "void rankoneinvit_runtime()" that tabulates
// the runtimes of both implementations according to the problem description.
// This subtask can be solved in many ways, and is not a part of the tests.
// In particular, the choice of the "tol" argument affects the runtime.
// Tip: To get your code to compile, start by defining all functions and
// include the appropriate header files. Then implement the functions.
// START

#include "timer.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

using namespace Eigen;

/* @brief Compute $l_{min}$ from vector $d$
 * Naive implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_0 */
double rankoneinvit(const VectorXd &d, const double &tol) {
  // Initialization
  double lmin = 0;
  VectorXd ev = d;
  double lnew = d.cwiseAbs().minCoeff();

  while (std::abs(lnew - lmin) > tol * lmin) {
    lmin = lnew;
    MatrixXd M = d.asDiagonal();
    M += ev * ev.transpose();
    // Instead of calculating the inverse of M directly
    // to get ev_new = M^{-1}*ev_old,
    // we solve the LSE M*ev_new = ev_old.
    ev = M.lu().solve(ev);
    // Equivalent to ev <- ev/|ev|:
    ev.normalize();
    lnew = ev.transpose() * M * ev;
  }

  return lnew;
}
/* SAM_LISTING_END_0 */

/* @brief Compute $l_{min}$ from vector $d$
 * Optimized implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_1 */
double rankoneinvit_fast(const VectorXd &d, const double &tol) {
  // Initialization
  double lmin = 0;
  VectorXd ev = d;
  double lnew = d.cwiseAbs().minCoeff();

  VectorXd dinv = (1 / d.array()).matrix();
  while (std::abs(lnew - lmin) > tol * lmin) {
    lmin = lnew;
    VectorXd ev0 = ev;

    // Here we solve the linear system
    // with the Sherman-Morrison-Woodbury formula
    // in the case of rank-1 perturbations.
    // This holds from $M = diag(d) + ev*ev^t$
    VectorXd Aib = dinv.cwiseProduct(ev);
    double temp = ev.transpose() * Aib;
    ev = Aib / (1 + temp);
    ev.normalize();
    // Better than the corresponding naive implementation.
    // This holds from $M = diag(d) + ev*ev^t$, too
    lnew = ev.transpose() * d.cwiseProduct(ev) + pow(ev.transpose() * ev0, 2);
  }

  return lnew;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void rankoneinvit_runtime() {
  unsigned int repeats = 3;
  double tol = 1e-3;
  double lmin;
  VectorXd d;
  Timer tm_slow, tm_fast;

  std::cout << std::endl
            << std::setw(15) << "n" << std::setw(15) << "Slow" << std::setw(15)
            << "Fast" << std::endl;

  for (unsigned int p = 2; p <= 8; p++) {
    tm_slow.reset();
    tm_fast.reset();
    unsigned int n = pow(2, p);

    for (unsigned int r = 0; r < repeats; ++r) {
      d = VectorXd::LinSpaced(n, 1, 2);

      tm_slow.start();
      lmin = rankoneinvit(d, tol);
      tm_slow.stop();

      tm_fast.start();
      lmin = rankoneinvit_fast(d, tol);
      tm_fast.stop();
    }
    std::cout << std::setw(15) << n << std::scientific << std::setprecision(3)
              << std::setw(15) << tm_slow.min() << std::setw(15)
              << tm_fast.min() << std::endl;
  }
}
/* SAM_LISTING_END_2 */

#endif
// END
