#ifndef TOEPLITZ_HPP
#define TOEPLITZ_HPP

#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "pconvfft.hpp"

using namespace Eigen;

/* @brief Build a Toeplitz matrix $\VT$ from $\Vc$ and $\Vr$
 * @param[in] c An $m$-dimensional vector, first column of $\VT$
 * @param[in] r An $n$-dimensional vector, first row of $\VT$
 * @param[out] T The $m \times n$ Toeplitz matrix from $\Vc$ and $\Vr$
 */
/* SAM_LISTING_BEGIN_5 */
MatrixXd toeplitz(const VectorXd & c, const VectorXd & r)
{
	if(c(0) != r(0)) {
		std::cerr << "First entries of c and r are different!" <<
		std::endl << "We assign the first entry of c to the diagonal"
		<< std::endl;
	}
	
	// Initialization
	int m = c.size();
	int n = r.size();
	MatrixXd T(m, n);
	
	// TODO: build Toeplitz matrix $\VT$
	// START
	
	// END
	
	return T;
}
/* SAM_LISTING_END_5 */

/* @brief Do something...
 * @param[in] c An $n$-dimensional vector
 * @param[in] r An $n$-dimensional vector
 * @param[in] x An $n$-dimensional vector
 * @param[out] y An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd toepmatmult(const VectorXd & c, const VectorXd & r,
					 const VectorXd & x)
{
	assert(c.size() == r.size() &&
		   c.size() == x.size() &&
		   "c, r, x have different lengths!");
	
	MatrixXd T = toeplitz(c,r);
	
	VectorXd y = T*x;
	
	return y;
}
/* SAM_LISTING_END_0 */

/* @brief Do something...
 * @param[in] c An $n$-dimensional vector
 * @param[in] r An $n$-dimensional vector
 * @param[in] x An $n$-dimensional vector
 * @param[out] y An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd toepmult(const VectorXd & c, const VectorXd & r,
				  const VectorXd & x)
{
	assert(c.size() == r.size() &&
		   c.size() == x.size() &&
		   "c, r, x have different lengths!");
	int n = c.size();
	
	VectorXcd cr_tmp = c.cast<std::complex<double>>();
	cr_tmp.conservativeResize(2*n); cr_tmp.tail(n) = VectorXcd::Zero(n);
	cr_tmp.tail(n-1).real() = r.tail(n-1).reverse();
	
	VectorXcd  x_tmp = x.cast<std::complex<double>>();
	x_tmp.conservativeResize(2*n);   x_tmp.tail(n) = VectorXcd::Zero(n);
	
	VectorXd y = pconvfft(cr_tmp, x_tmp).real();
	y.conservativeResize(n);
	
	return y;
}
/* SAM_LISTING_END_1 */

/* @brief Do something...
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd ttmatsolve(const VectorXd & h, const VectorXd & y)
{
	assert(h.size() == y.size() &&
		   "h and y have different lengths!");
	int n = h.size();
	
	VectorXd h_tmp = VectorXd::Zero(n);
	h_tmp(0) = h(0);
	
	MatrixXd T = toeplitz(h,h_tmp);
	
	VectorXd x = T.fullPivLu().solve(y);
	
	return x;
}
/* SAM_LISTING_END_2 */

/* @brief Do something...
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[in] l An integer
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_3 */
VectorXd ttrecsolve(const VectorXd & h, const VectorXd & y, int l)
{
	assert(h.size() == y.size() &&
		   "h and y have different lengths!");
	
	VectorXd x;
	
	if(l == 0) {
		x.resize(1);
		x(0) = y(0)/h(0);
	} else {
		int n = std::pow(2,l);
		int m = n/2;
		
		assert(h.size() == n && y.size() == n &&
			   "h and y have length different from 2^l!");
		
		VectorXd x1 = ttrecsolve(h.head(m), y.head(m), l-1);
		VectorXd y2 = y.segment(m,m) - toepmult(h.segment(m,m),
												h.segment(1,m).reverse(), x1);
		VectorXd x2 = ttrecsolve(h.head(m), y2, l-1);
		
		x.resize(n);
		x.head(m) = x1;
		x.tail(m) = x2;
	}
	
	return x;
}
/* SAM_LISTING_END_3 */

/* @brief Wrapper for 'ttrecsolve' for any size $n$
 * @param[in] h An $n$-dimensional vector
 * @param[in] y An $n$-dimensional vector
 * @param[out] x An $n$-dimensional vector
 */
/* SAM_LISTING_BEGIN_4 */
VectorXd ttsolve(const VectorXd & h, const VectorXd & y)
{
	assert(h.size() == y.size() &&
		   "h and y have different lengths!");
	int n = h.size();
	
	VectorXd x;
	
	// TODO: wrap 'ttrecsolve' for any size $n$
	// START
	
	// END
	
	return x;
}
/* SAM_LISTING_END_4 */

#endif
