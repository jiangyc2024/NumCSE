#include <cmath>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

VectorXcd fftreal(const VectorXd& y) {
	  
	size_t n = y.size();
	size_t m = n/2;
	assert(n % 2 != 0 &&
	"n must be even");
	
	VectorXcd y_(m);
	for(size_t j=0; j<m; ++j) {
		
		std::complex<double> tmp( y[2*j], y[2*j+1] );
		y_[i] = tmp;
	}
	
	FFT<double> fft;
	VectorXcd h;
	fft.fwd(h, y_);
	h.conservativeResize(m+1);
	h[m] = h[0];
	
	VectorXcd c(h.size());
	for(size_t j=0; j<(m+1); ++j) {
		
		c[j] = 0.5*(h[j] + h[m-j].conjugate());
		std::complex<double> tmp1( sin(2*M_PI*j/n), cos(2*M_PI*j/n) );
		std::complex<double> tmp2;
		tmp2 = 0.5*(h[j] - h[m-j].conjugate());
		c[j] = c[j] - tmp1*tmp2;
	}
	
	c.conservativeResize(2*m);
	for(size_t j=1; j<m; ++j) {
		c[m+j] = c[m-j];
	}

	return c;
}
