#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
void sinetransform(const Eigen::VectorXd &y, Eigen::VectorXd& s)
{
	int n = y.rows() + 1;
	std::complex<double> i(0,1);

	// Prepare sine terms
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n-1, 1, n-1);
	Eigen::VectorXd sinevals = x.unaryExpr([&](double z){ return imag(std::pow(std::exp(i*M_PI/(double)n), z)); });

	// Transform coefficients
	Eigen::VectorXd yt(n);
	yt(0) = 0;
	yt.tail(n-1) = sinevals.array() * (y + y.reverse()).array() + 0.5*(y-y.reverse()).array();

	// FFT
	Eigen::VectorXcd c;
	Eigen::FFT<double> fft;
	fft.fwd(c,yt);

	s.resize(n);
	s(0) = sinevals.dot(y);

	for (int k=2; k<=n-1; ++k)
	{
		int j = k-1; // Shift index to consider indices starting from 0
		if (k%2==0) 
			s(j) = -c(k/2).imag();
		else 
			s(j) = s(j-2) + c((k-1)/2).real();
	}
}
/* SAM_LISTING_END_0 */
