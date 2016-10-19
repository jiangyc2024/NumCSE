# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXd;
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
void freqfilter(const VectorXd& y, int k, VectorXd& low, VectorXd& high) {
  const unsigned m = y.size()/2;
  
  Eigen::FFT<double> fft;
  VectorXcd c = fft.fwd(y);
  
  VectorXcd clow = c;
  for (int j = -k; j <= +k; ++j) {
	  clow(m+j) = 0;
  }
  
  VectorXcd chigh = c - clow;
  
  low = fft.inv(clow).real();
  high = fft.inv(chigh).real();
}
/* SAM_LISTING_END_0 */
