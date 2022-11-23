#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
void icosinetransform(const Eigen::VectorXd& c, Eigen::VectorXd& y)
{
	size_t n = c.size();
	
	std::complex<double> i(0,1);
	Eigen::VectorXcd c_1(n);
	c_1(0) = sqrt(2)*c(0);
	for(size_t j=1; j<n; ++j) {
		c_1(j) = pow(exp(-i*M_PI/(2*(double)n)), j) * c(j);
	}
	
	Eigen::VectorXcd c_2(2*n);
	c_2.head(n) = c_1;
	c_2(n) = 0;
	c_2.tail(n-1) = c_1.tail(n-1).reverse().conjugate();
	
	// FFT
	Eigen::VectorXd z;
	Eigen::FFT<double> fft;
	fft.inv(z,c_2);
	
	// To obtain the same result of Matlab,
	// shift the inverse FFT result by 1.
	Eigen::VectorXd y_(2*n);
	y_.head(2*n-1) = z.tail(2*n-1);
	y_(2*n-1) = z(0);
	
	y = 2*y_.head(n);
}
/* SAM_LISTING_END_0 */
