#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <Eigen/Dense>
#include <math.h>

#include "fft.hpp"

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
bool testNormEqMatrix(unsigned int n, unsigned int m){

  // TODO: Test if the two definitions, given in (5.2.4) and 
  // (5.2.7), of $A^T*A$ are equal. 
  // START 
  MatrixXd A(n,m+1);
  // Initializing the vectors described in (5.2.4)
  VectorXd n_linear = VectorXd::LinSpaced(n,0,n-1);
  VectorXd m_linear = VectorXd::LinSpaced(m+1,0,m);
  // building A as described in (5.2.4)
  A = 2.0*M_PI/n* n_linear * m_linear.transpose();
  A = A.array().cos().matrix();
  // Using the knowledge gained in (5-2.b) to 
  // initialize $A^T*A$ directly as a diagonal matrix.
  MatrixXd ATA = MatrixXd::Zero(m+1, m+1);
  ATA.diagonal() = VectorXd::Constant(m+1, n/2.0);
  ATA(0,0) = n;
  // Comparing the two approaches to calculate $A^T*A$.
  if((ATA-A.transpose()*A).norm() < 1e-10 * m * n)
    return true;
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
 * \return The coefficients of the trigonometric polinomial.
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd find_c(const VectorXd& d, unsigned int m) {
  unsigned int n = d.size();

  // We will use a real to complex, discrete Fourier transform.
  FFT<double> fft;
  VectorXd rhs;
  // TODO: find the coefficients
  // START 
  // Computing DFT of d
  VectorXcd fourier = fft.fwd(d);
  // Gathering what is important for the rhs.
  rhs = fourier.real().head(m);
  // Solving normal equation by inverting diagonal matrix
  rhs /= n;
  rhs.tail(m-1) *= 2;
  // END
  return rhs;
}
/* SAM_LISTING_END_1 */

#endif
