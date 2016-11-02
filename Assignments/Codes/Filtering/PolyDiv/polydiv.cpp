#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "timer.h"

using namespace Eigen;

/* @brief 
 * @param[in] 
 * @param[out] 
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd polyMult_naive(const VectorXd & u, const VectorXd & v)
{
    // Initialization
    unsigned m = u.size();
    unsigned n = v.size();
    unsigned dim = max(m, n);
    
    VectorXd u_tmp = u;
    u_tmp.conservativeResize(dim);
    VectorXd v_tmp = v;
    v_tmp.conservativeResize(dim);
    
    VectorXd uv(m+n-1); // Degree is (m-1) + (n-1)

#if SOLUTION
    for(unsigned i=0; i<uv.size(); ++i) {
		unsigned fst = max(0, i - n - 1);
		unsigned lst = min(i, n - 1);
		for(unsigned j=fst; j<=lst; ++j) {
			uv(i) += u_tmp(j) * v_tmp(i-j);
		}
	}
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return uv;
}
/* SAM_LISTING_END_0 */

/* @brief 
 * @param[in] 
 * @param[out] 
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd polyMult_fast(const VectorXd & u, const VectorXd & v)
{
    // Initialization
    unsigned m = u.size();
    unsigned n = v.size();
    unsigned dim = max(m, n);
    
    VectorXd u_tmp = u;
    u_tmp.conservativeResize(u.size() + n - 1);
    VectorXd v_tmp = v;
    v_tmp.conservativeResize(v.size() + m - 1);

#if SOLUTION
	Eigen::FFT<double> fft;
	VectorXcd tmp = ( fft.fwd(u_tmp) ).cwiseProduct( fft.fwd(v_tmp) );
	VectorXd uv = fft.inv(tmp).real();
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	return uv;
}
/* SAM_LISTING_END_1 */

/* @brief 
 * @param[in] 
 * @param[out] 
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd polyDiv(const VectorXd & uv, const VectorXd & u)
{
    // Initialization
    unsigned mn = uv.size();
    unsigned m = u.size();
    
    VectorXd u_tmp = u;
    u_tmp.conservativeResize(mn);

#if SOLUTION
	Eigen::FFT<double> fft;
	VectorXcd tmp = ( fft.fwd(uv) ).cwiseQuotient( fft.fwd(u_tmp) );
	VectorXd v = fft.inv(tmp).real();
#else // TEMPLATE
    // TODO: 
#endif // TEMPLATE

	v.conservativeResize(mn - m + 1); // (mn-1) - (m-1) + 1
	return v;
}
/* SAM_LISTING_END_2 */

int main() {
	// Initialization
	unsigned m = 4;
    unsigned n = 3;
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
	unsigned repeats = 3;
	VectorXd uv;
	
    // Compute runtimes of different multiplicators
    std::cout << "Runtime comparison of naive v efficient multiplicator"
			  << std::endl;

    // Header
    std::cout << std::setw(20) << "n"
              << std::setw(20) << "time naive [s]"
              << std::setw(20) << "time fast [s]"
              << std::endl;

    // Loop over matrix size
    for(unsigned p = 2; p <= 9; ++p) {
        // Timers
        Timer tm_naive, tm_fast;
        unsigned n = pow(2,p);

        // Repeat test many times
        for(unsigned r = 0; r < repeats; ++r) {
            u = VectorXd::Random(n);
            v = VectorXd::Random(n);

            // Compute runtime of naive multiplicator
            tm_naive.start();
            uv = polyMult_naive(u, v);
            tm_naive.stop();
            // Compute runtime of efficient multiplicator
            tm_fast.start();
            uv = polyMult_fast(u, v);
            tm_fast.stop();
        }
        // Print runtimes
        std::cout << std::setw(20) << n
                  << std::scientific << std::setprecision(3)
                  << std::setw(20) << tm_naive.min()
                  << std::setw(20) << tm_fast.min()
                  << std::endl;
    }
}
