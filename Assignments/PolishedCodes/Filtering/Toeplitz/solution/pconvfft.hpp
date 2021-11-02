#ifndef PCONVFFT_HPP
#define PCONVFFT_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd pconvfft(const Eigen::VectorXcd &u,
                          const Eigen::VectorXcd &x) {
  Eigen::FFT<double> fft;
  Eigen::VectorXcd tmp = (fft.fwd(u)).cwiseProduct(fft.fwd(x));

  return fft.inv(tmp);
}
/* SAM_LISTING_END_0 */

#endif
