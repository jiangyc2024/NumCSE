#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
void cosinetransform(const Eigen::VectorXd& y, Eigen::VectorXd& c)
{
	int n = y.size();
	
	Eigen::VectorXd y_(2*n);
	y_.head(n) = y;
	y_.tail(n) = y.reverse();
	
	// FFT
	Eigen::VectorXcd z;
	Eigen::FFT<double> fft;
	fft.fwd(z,y_);
	
	std::complex<double> i(0,1);
	c.resize(n);
	c(0) = z(0).real()/(2*sqrt(2));
	for(size_t j=1; j<n; ++j) {
		c(j) = (0.5 * pow(exp(-i*M_PI/(2*(double)n)), j) * z(j)).real();
	}
}
/* SAM_LISTING_END_0 */
