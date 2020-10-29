#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <math.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

/*!
 * \brief testNormEqMatrix Create the matrix $A^TA$
 * in two different ways and make sure they are
 * approximately equal.
 * \param n number of different measurements
 * \param m degree of the trigonometric polynomial to
 * be fitted
 */
/* SAM_LISTING_BEGIN_0 */
bool testNormEqMatrix(unsigned int n, unsigned int m) {
  // TODO: (5-2.c) Test if the two definitions, given in (5.2.4) and
  // (5.2.7), of $A^T*A$ are equal.
  // START

  // END
  return false;
}
/* SAM_LISTING_0 */

/*!
 * \brief find_c Find best trigonometric polynomial
 * passing trough distances $\mathbf{d}$.
 * Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 * \param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * \param m Degree of the trigonometric polynomial $p$.
 * \return The coefficients of the trigonometric polynomial.
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd find_c(const VectorXd& d, unsigned int m) {
  unsigned int n = d.size();

  // We will use a real to complex, discrete Fourier transform.
  FFT<double> fft;
  VectorXd rhs;
  // TODO: (5-2.d) find the coefficients
  // START

  // END
  return rhs;
}
/* SAM_LISTING_END_1 */

#endif
