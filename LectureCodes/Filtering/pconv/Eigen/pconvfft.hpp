#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd pconvfft(const Eigen::VectorXcd &u,
                          const Eigen::VectorXcd &x) {
  Eigen::FFT<double> fft;
  return fft.inv(((fft.fwd(u)).cwiseProduct(fft.fwd(x))).eval());
}
/* SAM_LISTING_END_0 */
