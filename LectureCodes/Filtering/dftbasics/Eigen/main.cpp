/* Demonstraction code for course Numerical Methods  for CSE, ETH Zurich
   FFT in Eigen
   @author Ralf Hiptmair
   @date August 2020
*/

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <unsupported/Eigen/FFT>

/* SAM_LISTING_BEGIN_0 */
int main(int /*argc*/, char ** /*argv*/) {
  using Comp = std::complex<double>;
  const Eigen::VectorXcd::Index n = 5;
  Eigen::VectorXcd y(n);
  Eigen::VectorXcd c(n);
  Eigen::VectorXcd x(n);
  y << Comp(1, 0), Comp(2, 1), Comp(3, 2), Comp(4, 3), Comp(5, 4);
  Eigen::FFT<double> fft;  // DFT transform object
  c = fft.fwd(y);          // DTF of y, see \cref{def:DFT}
  x = fft.inv(c);          // inverse DFT of c, see \eqref{eq:invdft}

  std::cout << "y = " << y.transpose() << std::endl
            << "c = " << c.transpose() << std::endl
            << "x = " << x.transpose() << std::endl;
  return 0;
}
/* SAM_LISTING_END_0 */
