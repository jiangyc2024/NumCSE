#include <algorithm>
#include <cstdlib>
#include <limits>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

/* @brief Polynomial multiplication -- naive implementation
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[in] v Vector of coefficients of polynomial $v$
 * @param[out] uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd polyMult_naive(const VectorXd & u, const VectorXd & v)
{
	// Initialization
	int m = u.size();
	int n = v.size();
	int dim = std::max(m, n);
	
	VectorXd u_tmp = u;
	u_tmp.conservativeResizeLike(VectorXd::Zero(dim));
	VectorXd v_tmp = v;
	v_tmp.conservativeResizeLike(VectorXd::Zero(dim));
	
	VectorXd uv(m+n-1); // Degree is (m-1) + (n-1)
	
	// TODO: multiply polynomials $u$ and $v$ naively
	// START
	for(int i=0; i<uv.size(); ++i) {
		int fst = std::max(0, i - (dim - 1));
		int lst = std::min(i, dim - 1);
		for(int j=fst; j<=lst; ++j) {
			uv(i) += u_tmp(j) * v_tmp(i-j);
		}
	}
	// END
	
	return uv;
}
/* SAM_LISTING_END_0 */

/* @brief Polynomial multiplication -- efficient implementation
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[in] v Vector of coefficients of polynomial $v$
 * @param[out] uv Vector of coefficients of polynomial $uv = u*v$
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd polyMult_fast(const VectorXd & u, const VectorXd & v)
{
	// Initialization
	int m = u.size();
	int n = v.size();
	
	VectorXd u_tmp = u;
	u_tmp.conservativeResizeLike(VectorXd::Zero(u.size() + n - 1));
	VectorXd v_tmp = v;
	v_tmp.conservativeResizeLike(VectorXd::Zero(v.size() + m - 1));
	
	VectorXd uv;
	
	// TODO: multiply polynomials $u$ and $v$ efficiently
	// START
	Eigen::FFT<double> fft;
	VectorXcd u_tmp_ = fft.fwd(u_tmp);
	VectorXcd v_tmp_ = fft.fwd(v_tmp);
	VectorXcd tmp = u_tmp_.cwiseProduct(v_tmp_);
	uv = fft.inv(tmp).real();
	// END
	
	return uv;
}
/* SAM_LISTING_END_1 */

/* @brief Polynomial division -- efficient implementation
 * @param[in] uv Vector of coefficients of polynomial $uv$
 * @param[in] u Vector of coefficients of polynomial $u$
 * @param[out] uv Vector of coefficients of polynomial $v$
 * from $uv = u*v$ (division without remainder)
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd polyDiv(const VectorXd & uv, const VectorXd & u)
{
	// Initialization
	int mn = uv.size();
	int m = u.size();
	int dim = std::max(mn, m);
	
	VectorXd uv_tmp = uv;
	uv_tmp.conservativeResizeLike(VectorXd::Zero(dim));
	VectorXd u_tmp  = u;
	u_tmp.conservativeResizeLike(VectorXd::Zero(dim));
	
	/* SAM_LISTING_BEGIN_3 */
	// TODO: make sure that $uv$ can be divided by $u$
	// START
	// Machine precision
	double epsilon = std::numeric_limits<double>::epsilon();
	for(int i=dim-1; i>=0; --i) { // Right-most value of longest vector
		// is not necessarily numerically different than 0...
		if(std::abs(uv_tmp(i)) > epsilon*uv.norm()) {
			// No problem
			break;
		} else if(std::abs(u_tmp(i)) > epsilon*u.norm()) {
			std::exit(EXIT_FAILURE);
		}
	}
	// END
	/* SAM_LISTING_END_3 */
	
	VectorXd v;
	
	// TODO: divide polynomials $uv$ and $u$ efficiently (no remainder)
	// START
	Eigen::FFT<double> fft;
	VectorXcd uv_tmp_ = fft.fwd(uv_tmp);
	VectorXcd u_tmp_ = fft.fwd(u_tmp);
	VectorXcd tmp = uv_tmp_.cwiseQuotient(u_tmp_);
	v = fft.inv(tmp).real();
	
	v.conservativeResizeLike(VectorXd::Zero(mn - m + 1));
	// (mn-1) - (m-1) + 1
	// END
	return v;
}
/* SAM_LISTING_END_2 */
