#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <Eigen/Dense>
#include <cassert>

//! @brief Evaluate the Legendre polynomials and its derivatives in vector $x$
//! using the 3-term recursion formulae.
//! @param[in] evaluation points
//! @param[out] matrices $Lx$ and $DLx$.
/* SAM_LISTING_BEGIN_0 */
void legvals(const Eigen::VectorXd& x, Eigen::MatrixXd& Lx,
             Eigen::MatrixXd& DLx) {
  const unsigned int n = Lx.cols() - 1;
  const unsigned int N = x.size();

  // TODO: (8-10.c) Fill matrices Lx and DLx with values of the legendre
  // polynomials at x, i.e. $\{P_k(x_j)\}_{jk}$ and $\{P'_k(x_j)\}_{jk}$ START

  // END
}
/* SAM_LISTING_END_0 */

//! @brief Evaluate $P_n(x)$ for a scalar $x$ and integer $n$.
/* SAM_LISTING_BEGIN_1 */
double Pnx(const double x, const int n) {
  Eigen::VectorXd Px(n + 1);

  // TODO: (8-10.d, optional) Evaluate $P_n(x)$ for a scalar $x$ and integer
  // $n$. START

  // END
  return Px(n);
}
/* SAM_LISTING_END_1 */

//! @brief Find the Gauss points using the secant method without regula falsi.
/* SAM_LISTING_BEGIN_2 */
Eigen::MatrixXd gaussPts(const unsigned int n, const double rtol = 1e-10,
                         const double atol = 1e-12) {
  Eigen::MatrixXd zeros = Eigen::MatrixXd::Zero(n, n);

  // TODO: (8-10.d) Find the zeros of the legendre polynomials using the secant
  // method. START

  // END
  return zeros;
}
/* SAM_LISTING_END_2 */

//! @brief Find the Gauss points using the secant method with regula falsi.
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd gaussPts_regulaFalsi(const unsigned int n,
                                     const double rtol = 1e-10,
                                     const double atol = 1e-12) {
  Eigen::MatrixXd zeros = Eigen::MatrixXd::Zero(n, n);

  // TODO: (8-10.f) Find the zeros of the legendre polynomials using the secant
  // method with regula falsi. START

  // END
  return zeros;
}
/* SAM_LISTING_END_3 */

#endif
