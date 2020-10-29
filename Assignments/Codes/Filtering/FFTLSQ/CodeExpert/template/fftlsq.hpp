#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <math.h>

#include <Eigen/Dense>
#include <iomanip>
#include <unsupported/Eigen/FFT>
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/*!
 * \brief eval_p Given polynomial coefficients, return value of polynomial
 * at $n$ equidistant points.
 *
 * \param p Coefficient vector of trigonometrix polynomial.
 * \param n Number of equidistant points at which to evaluate.
 * \return Value of polynomial $p$ at $2\pi i / n$.
 */
VectorXd eval_p(const VectorXd& c, const unsigned int n) {
  // Degree of polynomial
  const unsigned int m = c.size();

  VectorXd ret(n);
  // Loop over all points
  for (unsigned int i = 0; i < n; ++i) {
    double r = 0;
    // Loop over all coefficients
    for (unsigned int j = 0; j < m; ++j) {
      r += c(j) * std::cos(2 * M_PI * i * j / n);
    }
    ret(i) = r;
  }
  return ret;
}


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


/*!
 * \brief using implementation of find_c compute coefficients c_k and p* for m=1,2,3 
 * using matplotlibcpp's plt() create a plot showing the ellipse
 * in the same plot insert the curves descriped by p* of degrees m=1,2,3
 */
/* SAM_LISTING_BEGIN_2 */
void fitEllipse(void) {
    
  unsigned int n = 10;
  unsigned int m = 3;
  
  // Test points
  constexpr unsigned int npoints = 10;
  VectorXd d(npoints);
  d << 0.987214, 1.03579, 0.997689, 0.917471, 1.00474, 0.92209, 1.03517,
      1.08863, 0.904992, 0.956089;
  

  plt::figure();
  // TODO: (5-2.e) Tabulate the coefficients of find_c for m = 1, 2, 3,
  // plot the ellipse and also the curves of the trigonimetric polynomials
  // START
  
  // END
  plt::savefig("cx_out/orbit.png");

}

/* SAM_LISTING_END_2 */

#endif
