# pragma once
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x) {
  Eigen::FFT<double> fft;
  VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
  return fft.inv(tmp);
}
/* SAM_LISTING_END_0 */
