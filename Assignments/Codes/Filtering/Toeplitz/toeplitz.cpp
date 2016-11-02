#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

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
		T.col(i).tail(m-i) = c.tail(m-i);
	}
	for(int i=0; i<m; ++i) {
		T.row(i).tail(n-i-1) = r.tail(n-i-1);
	} // Do not reassign the diagonal

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
	
	VectorXd y;

#if SOLUTION
    MatrixXd T = toeplitz(c,r);
	
    y = T*x;
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

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
	
	VectorXd y;

#if SOLUTION
	VectorXd cr_tmp = c;
	cr_tmp.conservativeResize(2*n);
	cr_tmp.tail(n-1) = r.tail(n-1).reverse();
	
	VectorXd x_tmp = x;
	x_tmp.conservativeResize(2*n);

    y = pconvfft(cr_tmp, x_tmp);
	y.conservativeResize(n);
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

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
	
	VectorXd x;

#if SOLUTION
	VectorXd h_tmp(n);
	h_tmp(0) = h(0);

    MatrixXd T = toeplitz(h,h_tmp);
    
    x = T.fullPivLu().solve(y);
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return y;
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
		VectorXd y2 = y.segment(m,n) - toepmult(h.segment(m,n-m),
			h.segment(1,m).reverse(), x1);
		VectorXd x2 = ttrecsolve();
	}
	
	

#if SOLUTION
	VectorXd h_tmp(n);
	h_tmp(0) = h(0);

    MatrixXd T = toeplitz(h,h_tmp);
    
    x = T.fullPivLu().solve(y);
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return x;
}
/* SAM_LISTING_END_3 */

int main() {
	// Initialization
	int m = 4;
	int n = 3;
	VectorXd u(m);
	VectorXd v(n);
	u << 1, 2, 3, 4;
	v << 10, 20, 30;

	// Compute with both functions
	std::cout << "Check that all functions are correct" << std::endl;

	VectorXd uv_1 = polyMult_naive(u, v);
	std::cout << "Naive multiplicator: "
			  << std::endl << uv_1 << std::endl;

	VectorXd uv_2 = polyMult_fast(u, v);
	std::cout << "Efficient multiplicator: "
			  << std::endl << uv_2 << std::endl;

	std::cout << "Error = " << (uv_1 - uv_2).norm() << std::endl;

	VectorXd v_new = polyDiv(uv_2, u);
	std::cout << "Error of efficient division = " << (v - v_new).norm()
			  << std::endl;

	// Initialization
	int repeats = 3;
	VectorXd uv;
#if INTERNAL
    // sizes   will contain the size of the vectors
    // timings will contain the runtimes in seconds
    std::vector<double> sizes, timings_naive, timings_effic;
#endif

	// Compute runtimes of different multiplicators
	std::cout << "Runtime comparison of naive v efficient multiplicator"
			  << std::endl;

	// Header
	std::cout << std::setw(20) << "n"
			  << std::setw(20) << "time naive [s]"
			  << std::setw(20) << "time fast [s]"
			  << std::endl;

	// Loop over vector size
	for(int p = 2; p <= 15; ++p) {
		// Timers
		Timer tm_naive, tm_effic;
		int n = pow(2,p);

		// Repeat test many times
		for(int r = 0; r < repeats; ++r) {
			u = VectorXd::Random(n);
			v = VectorXd::Random(n);

			// Compute runtime of naive multiplicator
			tm_naive.start();
			uv = polyMult_naive(u, v);
			tm_naive.stop();
			// Compute runtime of efficient multiplicator
			tm_effic.start();
			uv = polyMult_fast(u, v);
			tm_effic.stop();
		}
		
#if INTERNAL
        // Save results in a vector
        sizes.push_back(n); // Save vector sizes
        timings_naive.push_back(tm_naive.min());
        timings_effic.push_back(tm_effic.min());
#endif
		
		// Print runtimes
		std::cout << std::setw(20) << n
				  << std::scientific << std::setprecision(3)
				  << std::setw(20) << tm_naive.min()
				  << std::setw(20) << tm_effic.min()
				  << std::endl;
	}
	
#if INTERNAL
    mgl::Figure fig;
    fig.title("Comparison of timings of polynomial multiplication");
    fig.ranges(1, 35000, 1e-8, 1e3);
    fig.setlog(true, true); // set loglog scale
    fig.plot(sizes, timings_naive, " r+").label("naive");
    fig.plot(sizes, timings_effic, " b+").label("efficient");
    fig.fplot("1e-9*x^3", "k|").label("O(n^3)");
    fig.fplot("1e-8*x^2", "k-").label("O(n^2)");
    fig.fplot("1e-7*x", "k:").label("O(n)");
    fig.xlabel("Vector size (n)");
    fig.ylabel("Time [s]");
    fig.legend(0, 1);
    fig.save("polydiv_comparison.eps");
    fig.save("polydiv_comparison.png");
#endif
}
