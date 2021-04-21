#ifndef QUADU_HPP
#define QUADU_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

/*!
 * \brief quadU Function implementing weighted Gauss quadrature
 * \param f Integrand, function handle
 * \param n Number of nodes
 * \return Value of integral
 */
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
double quadU(const Function& f, const unsigned int n) {
  // Value of integral
  double I = 0;

  // TODO: (8-6.i) Integrate f using weighted Gauss quadrature
  // START
  for (unsigned int j = 0; j < n; ++j) {
    // Weight
    const double w =
        M_PI / (n + 1) * std::pow(std::sin((j + 1) * M_PI / (n + 1)), 2);

    // Node
    const double xi = std::cos((j + 1.) / (n + 1) * M_PI);

    I += w * f(xi);
  }
  // END
  return I;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testQuadU(unsigned int nmax = 20) {
  // "Exact" value of integral
  constexpr double exact = 0.483296828976607;
  // Integrand
  auto f = [](double t) { return 1. / (2. + std::exp(3. * t)); };

  plt::figure();

  // TODO: (8-6.j) Tabulate and plot the quadrature error.
  // START
  // memory allocation for plot
  std::vector<double> err(nmax);
  std::vector<int> num_pts(nmax);

  std::cout << std::setw(5) << "n" << std::setw(20) << "Error" << std::setw(20)
            << "q" << std::endl;

  for (unsigned int n = 1; n <= nmax; ++n) {
    // Compute error
    err[n - 1] = std::abs(quadU(f, n) - exact);
    num_pts[n - 1] = n;

    if (n > 1) {
      // Print table
      std::cout << std::setw(5) << n << std::setw(20) << err[n - 1]
                << std::setw(20) << err[n - 1] / err[n - 2] << std::endl;
    }
  }
  // Error plot rendered by matplotlibcpp
  plt::semilogy(num_pts, err, "^", {{"label", "error"}});
  plt::xlim(0.5, nmax + 0.5);
  plt::xticks(num_pts);
  plt::xlabel("No. of quadrature nodes");
  plt::ylabel("|Error|");
  plt::title("Convergence of quadU");
  plt::grid("True");
  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_2 */

#endif
