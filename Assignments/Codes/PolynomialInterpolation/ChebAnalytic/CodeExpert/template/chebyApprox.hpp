#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include "intpolyval.hpp"

using Eigen::VectorXcd;
using Eigen::VectorXd;
typedef std::complex<double> complex;

/*!
 * \brief Approximate the length of the ellipse
 * $\gamma_p(\theta):=cos(\theta-i \log(\rho),\forall
 * 0\leq\theta\leq2\pi,\rho>1$ by sampling for N equidistant arguments \theta
 * \param rho, the parameter of the ellipse curve
 * \param N, number of equidistant points for approximation
 * \return Value of the length of the ellipse
 */
/* SAM_LISTING_BEGIN_0 */
double lengthEllipse(double rho, unsigned int N) {
  // TODO: (7-5.a) Approximate the length of ellipse given $\rho$ and N
  // equidistant points (in $\theta$)
  double length = 0.0;
  // START

  // END
  return length;
}
/* SAM_LISTING_END_0 */

/*!
 * \brief Returns an approximation of
 * $\min_\rho \frac{\vert\gamma_\rho\vert \max_{z \in \gamma_p} \vert f(z)
 * \vert}{\pi d([-1,1], \gamma)}$ and the $\rho$ for which the minimum is
 * attained
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> bestBound() {
  // TODO: (7-5.d) Compute (7.5.2) with d from (7.5.3) and
  // $\vert\gamma_\rho\vert$ from (7-5.a)
  // START
  return {0, 0};
  // END
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
  // TODO: (7-5.f)
  std::vector<double> errInf, upperBound;
  // START

  // END
  return {errInf, upperBound};
}
/* SAM_LISTING_END_2 */
