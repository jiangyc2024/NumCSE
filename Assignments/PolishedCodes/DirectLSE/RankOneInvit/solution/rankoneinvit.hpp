#ifndef RANKONEINVIT_HPP
#define RANKONEINVIT_HPP

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "timer.h"

/**
 * @brief Compute $l_{min}$ from vector $d$. Naive implementation.
 *
 * @param d An $n$-dimensional vector
 * @param tol Scalar of type 'double', the tolerance
 * @return double lmin
 */
/* SAM_LISTING_BEGIN_0 */
double rankoneinvit(const Eigen::VectorXd& d, const double& tol) {
  // Initialization
  double lmin = 0, lnew = 0;
  // TODO: (3-10.a) Port the pseudocode in algorithm (3.10.1) to C++-code.
  // START
  Eigen::VectorXd ev = d;
  lnew = d.cwiseAbs().minCoeff();

  while (std::abs(lnew - lmin) > tol * lmin) {
    lmin = lnew;
    Eigen::MatrixXd M = d.asDiagonal();
    M += ev * ev.transpose();
    // Instead of calculating the inverse of M directly
    // to get ev_new = M^{-1}*ev_old,
    // we solve the LSE M*ev_new = ev_old.
    ev = M.lu().solve(ev);
    // Equivalent to ev <- ev/|ev|:
    ev.normalize();
    lnew = ev.dot(M * ev);
  }
  // END
  return lnew;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Compute $l_{min}$ from vector $d$. Optimized implementation.
 *
 * @param d An $n$-dimensional vector
 * @param tol Scalar of type 'double', the tolerance
 * @return double lmin
 */
/* SAM_LISTING_BEGIN_1 */
double rankoneinvit_fast(const Eigen::VectorXd& d, const double& tol) {
  // Initialization
  double lmin = 0, lnew = 0;
  // TODO: (3-10.c) Rewrite rankoneinvit() to have much better asymptotic
  // complexity.
  // START
  Eigen::VectorXd ev = d;
  lnew = d.cwiseAbs().minCoeff();

  Eigen::VectorXd dinv = (1. / d.array()).matrix();
  while (std::abs(lnew - lmin) > tol * lmin) {
    lmin = lnew;
    Eigen::VectorXd ev0 = ev;

    // Here we solve the linear system
    // with the Sherman-Morrison-Woodbury formula
    // in the case of rank-1 perturbations.
    // This holds from $M = diag(d) + ev*ev^t$
    Eigen::VectorXd Aib = dinv.cwiseProduct(ev);
    const double temp = ev.transpose() * Aib;
    ev = Aib / (1 + temp);
    ev.normalize();
    // Better than the corresponding naive implementation.
    // This holds from $M = diag(d) + ev*ev^t$, too
    lnew = ev.transpose() * d.cwiseProduct(ev) + pow(ev.transpose() * ev0, 2);
  }
  // END
  return lnew;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Tabulates the runtimes of the two different implementations for
 * different sizes of $d$.
 *
 */
/* SAM_LISTING_BEGIN_2 */
void rankoneinvit_runtime() {
  constexpr unsigned int repeats = 3;
  constexpr double tol = 1e-3;
  // TODO: (3-10.e) Tabulate the runtimes of both implementations according to
  // the problem description. The choice of the "tol" argument affects the
  // runtime.
  // START
  double lmin;
  Eigen::VectorXd d;
  Timer tm_slow, tm_fast;

  std::cout << std::endl
            << std::setw(15) << "n" << std::setw(15) << "Slow" << std::setw(15)
            << "Fast" << std::endl;

  for (unsigned int n = 2; n <= 256; n <<= 1) {
    tm_slow.reset();
    tm_fast.reset();
    for (unsigned int r = 0; r < repeats; ++r) {
      d = Eigen::VectorXd::LinSpaced(n, 1, 2);

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
  // END
}
/* SAM_LISTING_END_2 */

#endif