#ifndef FFTLSQ_HPP
#define FFTLSQ_HPP

#include <Eigen/Dense>

#include "fft.hpp"

using namespace Eigen;

/*!
 * \brief gauss_fit Find best trigonometric polynomial
 * passing trough distances $\mathbf{d}$.
 * Using Least Squares formulation and assuming best fitting
 * polynomial has degree $m$.
 * \param d Vector of size $n$ distances at angle $2*\pi*i/n$.
 * \param m Degree of the trigonometric polynomial $p$.
 * \return The coefficients of the trigonometric polinomial.
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd gauss_fit(const VectorXd & d,
                   unsigned int m) {
    unsigned int n = d.size();

    // We will use a real to complex, discrete Fourier transform.
    FFT<double> fft;
	VectorXd rhs;
	// TODO: find the coefficients
	// START
    
	// END

    return rhs;
}
/* SAM_LISTING_END_1 */

#endif
