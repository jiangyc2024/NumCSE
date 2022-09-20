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
  // Method 1: $s = \sqrt(1 \pm t)$
  // Will use same node twice
  QuadRule Q = gauleg(n);

  // quadrature done on interval [0, 1], i.e. nodes and weights need to be
  // shifted
  for (unsigned i = 0; i < n; ++i) {
    // Transform nodes
    double x = (Q.nodes(i) + 1.) / 2.;
    // Weights
    double w = Q.weights(i) * x * x * std::sqrt(2. - x * x);
    // Symmetric summation
    I += w * (f(x * x - 1) + f(-x * x + 1));
  }

  // Method 2: $t = cos(s)$
  // We are actually using twice the number of nodes than Method 1
  //  QuadRule Q = gauleg(2 * n);
  //
  // quadrature done on interval [-pi/2, pi/2], i.e. nodes and weights need to
  // be shifted
  //  for(unsigned i = 0; i < 2 * n; ++i) {
  // Evaluate transformation
  //    double x = sin(Q.nodes(i) * M_PI_2);
  // Weights
  //    double w = Q.weights(i) * cos(Q.nodes(i) * M_PI_2) * cos(Q.nodes(i) *
  //    M_PI_2); I += w * f(x) * M_PI_2;
  //  }
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
  std::vector<int> num_nodes(max_N);  // number of nodes on x axis
  std::vector<double> err(max_N);     // error on y axis

  // Table header
  std::cout << std::setw(3) << "N" << std::setw(15) << "I_approx"
            << std::setw(15) << "error" << std::endl;

  for (unsigned N = 1; N <= max_N; ++N) {
    const double I_approx = quadsingint(f, N);
    num_nodes[N - 1] = N;
    err[N - 1] = std::abs(I_ex - I_approx);

    std::cout << std::setw(3)
              << N
              // Value of integral
              << std::setw(15)
              << I_approx
              // Error
              << std::setw(15) << err[N - 1] << std::endl;
  }
  plt::semilogy(num_nodes, err, "+", {{"label", "Error"}});
  plt::xlabel("No. of quadrature nodes");
  plt::ylabel("|Error|");
  plt::title("Convergence of quadsingint");
  plt::grid("True");
  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_2 */

#endif
