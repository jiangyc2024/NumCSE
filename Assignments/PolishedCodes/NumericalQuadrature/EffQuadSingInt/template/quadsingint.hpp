#ifndef QUADSINGINT_HPP
#define QUADSINGINT_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "gauleg.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//! @brief Compute integral $\int_{-1}^1 \sqrt(1-t^2) f(t) dt$ using
//! transformation
//! @param[in] n number of Gauss nodes (evaluate f at 2*n points)
//! @return value of integral
/* SAM_LISTING_BEGIN_1 */
template <class Function>
double quadsingint(const Function& f, const unsigned n) {
  double I = 0.;

  // TODO: (7-9.c) Compute $\int_{-1}^1 \sqrt(1-t^2) f(t) dt$ using
  // Gauss-Legendre quadrature. Ensure asymptotic exponential convergence. START

  // END
  return I;
}
/* SAM_LISTING_END_1 */

//! @brief Tabulates and plots the quadrature error, i.e. $\abs{W(f) -
//! \text{quadsingint(f,n)}}$ for $f(t) = \frac{1}{2 + e^{3t}}$ and $n =
//! 1,2,...,25$.
/* SAM_LISTING_BEGIN_2 */
void tabAndPlotQuadErr() {
  // Max num of Gauss points to use (in each direction)
  constexpr unsigned max_N = 25;
  // "Exact" integral
  constexpr double I_ex = 0.483296828976607;

  auto f = [](double t) { return 1. / (2. + std::exp(3. * t)); };

  plt::figure();
  // TODO: (7-9.e) Tabulate and plot quadrature error for the given function.
  // START

  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_2 */

#endif
