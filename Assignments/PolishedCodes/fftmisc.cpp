/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 * HW problem "Discrete Fourier Transform and Applications", prb:fftm
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <vector>

namespace FFTMisc {
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXcd inverseDFT(const Eigen::VectorXcd &c) {
  const Eigen::VectorXcd::Index n = c.size();
  Eigen::FFT<double> fft;
  Eigen::VectorXcd tmp = c.conjugate();
  Eigen::VectorXcd y = fft.fwd(tmp);
  return y.conjugate() / n;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXcd circmatvec(const Eigen::VectorXcd &u,
                            const Eigen::VectorXcd &x) {
  Eigen::FFT<double> fft;
  Eigen::VectorXcd tmp = (fft.fwd(u)).cwiseProduct(fft.fwd(x));
  return fft.inv(tmp);
}
/* SAM_LISTING_END_2 */

}  // namespace FFTMisc

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Inverse DFT via forward DFT" << std::endl;
  Eigen::VectorXcd c = Eigen::VectorXcd::Random(16);
  Eigen::VectorXcd y = FFTMisc::inverseDFT(c);
  Eigen::FFT<double> fft;
  Eigen::VectorXcd r = fft.fwd(y);
  std::cout << "Deviation = " << (c - r).norm() << std::endl;
  return 0;
}
