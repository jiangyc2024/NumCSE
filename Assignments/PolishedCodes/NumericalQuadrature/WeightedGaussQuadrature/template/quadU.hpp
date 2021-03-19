#ifndef QUADU_HPP
#define QUADU_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
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
template<typename Function>
double quadU(const Function& f, const unsigned int n) {
  // Value of integral
  double I = 0;
  
  // TODO: (8-6.i) Integrate f using weighted Gauss quadrature
  // START
  
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
  
  // END
  plt::savefig("./cx_out/convergence.png");
}
/* SAM_LISTING_END_2 */

#endif
