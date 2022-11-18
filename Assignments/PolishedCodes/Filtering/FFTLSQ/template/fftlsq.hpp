#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <math.h>

#include <Eigen/Dense>
#include <iomanip>
#include <unsupported/Eigen/FFT>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * @brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * @param c Coefficient vector of trigonometrix polynomial.
 * @param n Number of equidistant points at which to evaluate.
 * @return Eigen::VectorXd Value of polynomial $p$ at $2\pi i / n$.
 */
Eigen::VectorXd eval_p(const Eigen::VectorXd &c, const unsigned int n) {
  // Degree of polynomial
  const unsigned int m = c.size();

  Eigen::VectorXd ret(n);
  // Loop over all points
  for (unsigned int j = 0; j < n; ++j) {
    double r = 0;
    // Loop over all coefficients
    for (unsigned int k = 0; k < m; ++k) {
      r += c(k) * std::cos(2.0 * M_PI * k * j / n);
    }
    ret(j) = r;
  }

  return ret;
}

/**
 * @brief testNormEqMatrix Create the matrix $A^TA$ in two different ways and
 * make sure they are approximately equal.
 *
 * @param n number of different measurements
 * @param m degree of the trigonometric polynomial to
 * be fitted
 * @return bool indicating if test passed or failed.
 */
/* SAM_LISTING_BEGIN_0 */
bool testNormEqMatrix(unsigned int n, unsigned int m) {
  // TODO: (4-2.c) Test if the two definitions, given in (4.2.4) and
  // (4.2.7), of $A^T*A$ are equal.

  // START

  // END

  return false;
}
/* SAM_LISTING_END_0 */

/**
 * @brief find_c Find Best trigonometric polynomial passing trough distances
 * $\mathbf{d}$. Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 *
 * @param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * @param m Degree of the trigonometric polynomial $p$.
 * @return Eigen::VectorXd The coefficients of the trigonometric polynomial.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd find_c(const Eigen::VectorXd &d, unsigned int m) {
  unsigned int n = d.size();

  // We will use a real to complex, discrete Fourier transform.
  Eigen::FFT<double> fft;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m + 1);

  // TODO: (4-2.d) find the coefficients
  // START

  // END

  return rhs;
}
/* SAM_LISTING_END_1 */

/**
 * @brief using implementation of find_c compute coefficients c_k and p* for
 * m=1,2,3 using matplotlibcpp's plt() create a plot showing the ellipse in the
 * same plot insert the curves descriped by p* of degrees m=1,2,3
 *
 */
/* SAM_LISTING_BEGIN_2 */
void fitEllipse() {
  unsigned int n = 10;
  unsigned int m = 3;

  // Test points
  constexpr unsigned int npoints = 10;
  Eigen::VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;

  plt::figure();

  // TODO: (4-2.e) Tabulate the coefficients of find\_c for all $m=1,2,3$
  // plot the ellipse and also the curves of the trigonimetric polynomials
  // START

  // END
  plt::savefig("cx_out/orbit.png");
}
/* SAM_LISTING_END_2 */

#endif
