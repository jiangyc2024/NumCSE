#ifndef CHEBYAPPROX_HPP
#define CHEBYAPPROX_HPP

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "intpolyval.hpp"

using complex = std::complex<double>;

/**
 * @brief Approximate the length of the ellipse
 * $\gamma_p(\theta):=cos(\theta-i \log(\rho),\forall
 * 0\leq\theta\leq2\pi,\rho>1$ by sampling for N equidistant arguments \theta
 *
 * @param rho the parameter of the ellipse curve
 * @param N number of equidistant points for approximation
 * @return double Value of the length of the ellipse
 */
/* SAM_LISTING_BEGIN_0 */
double lengthEllipse(double rho, unsigned int N) {
  double length = 0.0;
  // TODO: (6-5.a) Approximate the length of ellipse given $\rho$ and N
  // equidistant points (in $\theta$)
  // START

  // END
  return length;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Returns an approximation of
 * $\min_\rho \frac{\vert\gamma_\rho\vert \max_{z \in \gamma_p} \vert f(z)
 * \vert}{\pi d([-1,1], \gamma)}$ and the $\rho$ for which the minimum is
 * attained
 *
 * @return std::pair<double, double> pair (min, argmin)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> bestBound() {
  const double maxRho = exp(asinh(M_PI / 3.));
  constexpr unsigned int N = 1000;
  double m = 0.0, r = 0.0;  // min and rho for output
  // TODO: (6-5.d) Compute (6.5.4) with d from (6.5.5) and
  // $\vert\gamma_\rho\vert$ from (6-5.a)
  // START

  // END
  return std::make_pair(m, r);
}
/* SAM_LISTING_END_1 */

/*!
 * \brief Tabulate the estimated upperbound and actual error using 10000 points
 * for LInf error of Chebychev interpolation with order from 4 to n_max
 * \param n_max the maximum order for interpolation in the table
 * \return the pair of the list of calculated error and the list of estimated
 * upperbound
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<std::vector<double>, std::vector<double>> compareErrorAndBound(
    unsigned int n_max) {
  std::vector<double> upperBound, errInf;
  constexpr unsigned N = 10000;

  // TODO: (6-5.f) Tabulate the L-inf error and the bounds.
  // START

  // END
  return {errInf, upperBound};
}
/* SAM_LISTING_END_2 */

#endif