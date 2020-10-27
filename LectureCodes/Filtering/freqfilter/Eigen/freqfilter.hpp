#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
void freqfilter(const VectorXd &y, int k, VectorXd &low, VectorXd &high) {
  const VectorXd::Index n = y.size();
  if (n % 2 != 0)
    throw std::runtime_error("Even vector length required!");
  const VectorXd::Index m = y.size() / 2;

  Eigen::FFT<double> fft;   // DFT helper object
  VectorXcd c = fft.fwd(y); // Perform DFT of input vector

  VectorXcd clow = c;
  // Set high frequency coefficients to zero, \cref{speccircle}
  for (int j = -k; j <= +k; ++j) clow(m + j) = 0;
  // (Complementary) vector of high frequency coefficients
  VectorXcd chigh = c - clow;

  // Recover filtered time-domain signals
  low = fft.inv(clow).real();
  high = fft.inv(chigh).real();
}
/* SAM_LISTING_END_0 */
