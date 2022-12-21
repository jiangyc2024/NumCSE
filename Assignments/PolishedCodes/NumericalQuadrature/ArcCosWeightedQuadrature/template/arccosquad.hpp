#ifndef ARCCOSQUAD_HPP
#define ARCCOSQUAD_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// The function gaussquad() for pre-computed (and pre-compiled)
// Gauss-Legendre quadrature rules up to 256 nodes.
#include "gaussquad.hpp"

/* SAM_LISTING_BEGIN_0 */
struct cvgdata {
  long n_{0};          // number of quadrature points
  double val_{0.0};    // approximate value
  double err_{-1.0};   // quadrature error
  double rate_{-1.0};  // estimate for rate of algebraic convergence
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void testConvGaussQuad() {
  // Table header
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10) << (" rate\n " + std::string(65, '-'))
            << std::endl;

  // Value obtained by overkill quadrature as reference value
  constexpr double ref_val = 1.7576137811123187;

  // TODO: (7-13.a) Print a table that allows you to predict the asymptotic
  // behaviour of Gauss-Legendre numerical quadrature when approximating I(f).
  // START

  // END
}
/* SAM_LISTING_END_1 */

/**
 * \brief Calculates the quadrature of f using transformation.
 *
 * \tparam FUNCTION suitable functor
 * \param f function to be integrated on [-1, 1]
 * \param n number of evaluations
 * \return double integral
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTION>
double arccosWeightedQuad(FUNCTION&& f, unsigned int n) {
  double I = 0.0;  // For accumulating quadrature result
  // TODO: (7-13.c) Approximate I(f) with exponential convergence in n.
  // START

  // END
  return I;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void testConvTrfGaussQuad() {
  // Table header
  std::cout << std::setw(5) << "n" << std::setw(20) << " I_n " << std::setw(30)
            << " Error " << std::setw(10)
            << (" decay ratio\n " + std::string(65, '-')) << std::endl;

  // Value obtained by overkill quadrature as reference value
  constexpr double ref_val = 1.7576137811123187;

  // TODO: (7-13.d) Print a table that allows you to predict the asymptotic
  // behaviour of arccosWeightedQuad().
  // START

  // END
}
/* SAM_LISTING_END_3 */

#endif