#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <unsupported/Eigen/FFT>

namespace fftreal {


using Eigen::VectorXd;
using Eigen::VectorXcd;

inline
/* SAM_LISTING_BEGIN_0 */
//  Perform fft on a \com{real} vector \Blue{$\Vy$} of even
//  length and return (complex) coefficients in \Blue{$\Vc$}
//  Note: \eigen's DFT method fwd() has this already implemented and
//        we could also just call: c = fft.fwd(y);
void fftreal(const VectorXd& y, VectorXcd& c) {
  const Eigen::Index n = y.size();
  const Eigen::Index m = n/2;
  if (n % 2 != 0) { std::cout << "n must be even!\n"; return; }
  
  // Step I: compute \Blue{$\Vh$} from \eqref{rft:1}, \eqref{rft:2}
  const std::complex<double> i(0,1); // Imaginary unit
  VectorXcd yc(m);
  for (Eigen::Index j = 0; j < m; ++j) {
    yc(j) = y(2*j) + i*y(2*j + 1);
  }
  
  Eigen::FFT<double> fft;
  VectorXcd d = fft.fwd(yc);
  VectorXcd h(m + 1);
  h << d, d(0);
  
  c.resize(n);
  // Step II: implementation of \eqref{rft:4}
  for (Eigen::Index k = 0; k < m; ++k) {
    c(k) = (h(k) + std::conj(h(m-k)))/2. - i/2.*std::exp(-2.*static_cast<double>(k)/static_cast<double>(n)*M_PI*i)*(h(k) - std::conj(h(m-k)));
  }
  c(m) = std::real(h(0)) - std::imag(h(0));
  for (Eigen::Index k = m+1; k < n; ++k) {
    c(k) = std::conj(c(n-k));
  }
}
/* SAM_LISTING_END_0 */


} //namespace fftreal