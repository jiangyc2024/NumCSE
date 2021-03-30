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
  assert(n >= 2);
  for (unsigned int j = 0; j < N; ++j) {
    Lx(j, 0) = 1.;
    Lx(j, 1) = x(j);
    DLx(j, 0) = 0;
    DLx(j, 1) = 1.;
    for (unsigned int k = 2; k <= n; ++k) {
      Lx(j, k) =
          (2 * k - 1.) / k * x(j) * Lx(j, k - 1) - (k - 1.) / k * Lx(j, k - 2);
      DLx(j, k) = (2 * k - 1.) / k * Lx(j, k - 1) +
                  (2 * k - 1.) / k * x(j) * DLx(j, k - 1) -
                  (k - 1.) / k * DLx(j, k - 2);
    }
  }
  // END
}
/* SAM_LISTING_END_0 */

//! @brief Evaluate $P_n(x)$ for a scalar $x$ and integer $n$.
/* SAM_LISTING_BEGIN_1 */
double Pnx(const double x, const int n) {
  Eigen::VectorXd Px(n + 1);

  // TODO: (8-10.d, optional) Evaluate $P_n(x)$ for a scalar $x$ and integer
  // $n$. START
  assert(n >= 1);
  Px(0) = 1.;
  Px(1) = x;
  for (unsigned int k = 2; k <= n; ++k) {
    Px(k) = (2 * k - 1.) / k * x * Px(k - 1) - (k - 1.) / k * Px(k - 2);
  }
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
  double x0, x1, f0, f1, s;
  for (unsigned int k = 1; k <= n; ++k) {
    for (unsigned int j = 1; j <= k; ++j) {
      // Initial guesses
      if (j == 1)
        x0 = -1.;
      else
        x0 = zeros(j - 2, k - 2);
      if (j == k)
        x1 = 1.;
      else
        x1 = zeros(j - 1, k - 2);

      // Secant method
      f0 = Pnx(x0, k);
      for (unsigned int i = 0; i < 1e4; ++i) {
        f1 = Pnx(x1, k);
        s = f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        f0 = f1;
        x1 -= s;
        if ((std::abs(s) <
             std::max(atol, rtol * std::min(std::abs(x0), std::abs(x1))))) {
          zeros(j - 1, k - 1) = x1;
          break;
        }
      }
    }
  }
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
  double x0, x1, f0, f1, s;
  for (unsigned int k = 1; k <= n; ++k) {
    for (unsigned int j = 1; j <= k; ++j) {
      // Initial guesses
      if (j == 1)
        x0 = -1;
      else
        x0 = zeros(j - 2, k - 2);
      if (j == k)
        x1 = 1;
      else
        x1 = zeros(j - 1, k - 2);

      // Secant method
      f0 = Pnx(x0, k);
      for (unsigned int i = 0; i < 1e4; ++i) {
        f1 = Pnx(x1, k);
        s = f1 * (x1 - x0) / (f1 - f0);
        if (Pnx(x1 - s, k) * f1 < 0) {  // regula falsi
          x0 = x1;
          f0 = f1;
        }
        x1 -= s;
        if ((std::abs(s) <
             std::max(atol, rtol * std::min(std::abs(x0), std::abs(x1))))) {
          zeros(j - 1, k - 1) = x1;
          break;
        }
      }
    }
  }
  // END
  return zeros;
}
/* SAM_LISTING_END_3 */

#endif
