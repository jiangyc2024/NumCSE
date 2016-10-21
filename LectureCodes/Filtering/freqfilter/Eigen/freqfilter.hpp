# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXd;
using Eigen::VectorXcd;

void freqfilter(const VectorXd& y, const int k, VectorXd& low, VectorXd& high) {
  const int m = y.size()/2;
  Eigen::FFT<double> fft;
  VectorXcd c = fft.fwd(y);
  VectorXcd clow = c; 
  clow.segment(m-k,m+k).setZero();
  VectorXcd chigh = c - clow;
  low = fft.inv(clow).real();
  high = fft.inv(chigh).real();
}
