#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
void sinetransformwraparound(const Eigen::VectorXd &y, Eigen::VectorXd& c)
{
	int n = y.rows()+1;
	Eigen::VectorXd yt(2*n);
	yt << 0, y, 0, -y.reverse();

	Eigen::VectorXcd ct;
	Eigen::FFT<double> fft;
	fft.SetFlag(Eigen::FFT<double>::Flag::Unscaled);
	fft.fwd(ct,yt);

	std::complex<double> v(0,2);
	c = (-ct.middleRows(1, n-1) / v).real();
}
/* SAM_LISTING_END_0 */
