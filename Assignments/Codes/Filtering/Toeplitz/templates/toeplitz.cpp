#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "timer.h"

#include "pconvfft.hpp"

using namespace Eigen;

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
 */
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
    
	for(int i=0; i<n; ++i) {
		T.col(i).tail(m-i) = c.head(m-i);
	}
	for(int i=0; i<m; ++i) {
		T.row(i).tail(n-i-1) = r.segment(1,n-i-1);
	} // Do not reassign the diagonal!

	return T;
}

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
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

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
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
	cr_tmp.conservativeResize(2*n);
	cr_tmp.tail(n-1).real() = r.tail(n-1).reverse();
	
	VectorXcd  x_tmp = x.cast<std::complex<double>>();
	x_tmp.conservativeResize(2*n);

    VectorXd y = pconvfft(cr_tmp, x_tmp).real();
	y.conservativeResize(n);

	return y;
}
/* SAM_LISTING_END_1 */

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd ttmatsolve(const VectorXd & h, const VectorXd & y)
{
	assert(h.size() == y.size() &&
		   "h and y have different lengths!");
	int n = h.size();
	
	VectorXd h_tmp(n);
	h_tmp(0) = h(0);

    MatrixXd T = toeplitz(h,h_tmp);
    
    VectorXd x = T.fullPivLu().solve(y);

	return x;
}
/* SAM_LISTING_END_2 */

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
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

/* @brief 
 * @param[in] 
 * @param[in] 
 * @param[out] 
 */
/* SAM_LISTING_BEGIN_4 */
VectorXd ttsolve(const VectorXd & h, const VectorXd & y)
{
	assert(h.size() == y.size() &&
		   "h and y have different lengths!");
	int n = h.size();
		   
	VectorXd x;
	
    // TODO: 

	return x;
}
/* SAM_LISTING_END_4 */

int main() {
	int n;
	
	// Initialization
	n = 3;
	VectorXd c(n), r(n), x(n);
	c << 1, 2, 3;
	r << 1, 4, 5;
	x << 6, 7, 8;

	// Compute with both functions toepmatmult and toepmult
	std::cout << "Check that toepmatmult and toepmult are correct"
			  << std::endl;
	VectorXd y_1 = toepmatmult(c, r, x);
	VectorXd y_2 = toepmult(c, r, x);
	std::cout << "Error = " << (y_1 - y_2).norm() << std::endl;

	// Initialization
	n = 4;
	VectorXd h(n), y(n);
	h << 1, 2, 3, 4;
	y << 5, 6, 7, 8;

	// Compute with both functions ttmatsolve and ttrecsolve
	std::cout << "Check that ttmatsolve and ttrecsolve are correct"
			  << std::endl;
	VectorXd x_1 = ttmatsolve(h, y);
	VectorXd x_2 = ttrecsolve(h, y, 2);
	std::cout << "Error = " << (x_1 - x_2).norm() << std::endl;

	// Initialization
	int repeats = 3;
	VectorXd out;

	// Compute runtimes of different multiplicators
	std::cout << "Runtime comparison of "
			  << "toepmatmult vs toepmult and ttmatsolve vs ttrecsolve"
			  << std::endl;

	// Loop over vector size
	for(int l = 3; l <= 11; ++l) {
		// Timers
		Timer tm_matmult, tm_mult, tm_ttmat, tm_ttrec;
		int n = pow(2,l);

		// Repeat test many times
		for(int repeat = 0; repeat < repeats; ++repeat) {
			c = VectorXd::Random(n);
			r = VectorXd::Random(n); r(0) = c(0);
			x = VectorXd::Random(n);
			h = VectorXd::LinSpaced(n,1,n).cwiseInverse();
			y = VectorXd::Random(n);

			// Compute runtime of toepmatmult
			tm_matmult.start();
			out = toepmatmult(c, r, x);
			tm_matmult.stop();
			// Compute runtime of toepmult
			tm_mult.start();
			out = toepmult(c, r, x);
			tm_mult.stop();
			// Compute runtime of ttmatsolve
			tm_ttmat.start();
			out = ttmatsolve(h, y);
			tm_ttmat.stop();
			// Compute runtime of ttrecsolve
			tm_ttrec.start();
			out = ttrecsolve(h, y, l);
			tm_ttrec.stop();
		}
		
		
		// Print progress
		std::cout << n << " completed" << std::endl;
	}
	
}
