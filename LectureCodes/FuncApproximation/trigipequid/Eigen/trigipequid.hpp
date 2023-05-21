# include <Eigen/Dense>
# include <cmath> 
# include <complex>
# include <iostream>
# include <unsupported/Eigen/FFT>

namespace trigipequid {


using Eigen::VectorXd;
using Eigen::VectorXcd; // complex Eigen vector

inline
/* SAM_LISTING_BEGIN_0 */
// Efficient FFT-based computation of coefficients in expansion \eqref{eq:trigpreal}
// for a trigonometric interpolation polynomial in equidistant points
// \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}.
// IN : \texttt{y} has to be a row vector of odd length, return values are column vectors
// OUT: vectors \Blue{$\left[\alpha_j\right]_j$}, \Blue{$\left[\alpha_j\right]_j$} of expansion coefficients
// with respect to trigonometric basis from \cref{def:trip}
std::pair<VectorXd,VectorXd> trigipequid(const VectorXd& y) {
  using index_t = VectorXcd::Index;
  const index_t N = y.size();
  const index_t n = (N - 1)/2;
  if (N % 2 != 1) {
    throw std::logic_error("Number of points must be odd!");
  }
  // prepare data for fft
  const std::complex<double> M_I(0,1); // imaginary unit
  // right hand side vector \Blue{$\Vb$} from \eqref{tip:FM}
  VectorXcd b(N);
  for (index_t k = 0; k < N; ++k) {
    b(k) = y(k) * std::exp(2*M_PI*M_I*(static_cast<double>(n)/static_cast<double>(N*k))); 
  }
  Eigen::FFT<double> fft; // DFT helper class
  VectorXcd c = fft.fwd(b); // means that ``c = fft(b)''

  // By \eqref{eq:trigpcomp} we can recover
  // \Blue{$\alpha_j = \frac{1}{2}(\gamma_{n-j}+\gamma_{n+j})$} and \Blue{$\beta_j = \frac{1}{2i}(\gamma_{n-j}-\gamma_{n+j})$}, \Blue{$j=1,\ldots,n$}, \Blue{$\alpha_0 = \gamma_n$}.
  VectorXcd alpha(n + 1);
  VectorXcd beta(n);
  alpha(0) = c(n)/static_cast<double>(N);
  for (index_t l = 1; l <= n; ++l) {
    alpha(l) = (c(n-l)+c(n+l))/static_cast<double>(N);
    beta(l-1) = -M_I*(c(n-l)-c(n+l))/static_cast<double>(N);
  }
  return {alpha.real(),beta.real()};
}
/* SAM_LISTING_END_0 */

inline
/* SAM_LISTING_BEGIN_1 */
// Efficient FFT-based computation of coefficients in expansion \eqref{eq:trigpreal}
// for a trigonometric interpolation polynomial in equidistant points
// \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}.
// IN : \texttt{y} has to be a row vector of odd length, return values are column vectors
//      \texttt{a}, \texttt{b} will be used to save the expansion coefficients
void trigipequid(const VectorXd& y, VectorXcd& a, VectorXcd& b) {
  const unsigned N = y.size();
  if (N % 2 != 1) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const unsigned n = (N - 1)/2;
  // prepare data for fft
  const std::complex< double > i(0,1); // imaginary unit
  std::vector< std::complex<double> > f(N);
  std::vector< std::complex<double> > c;

  for (unsigned k = 0; k < N; ++k) {
    // see \eqref{tip:FM}
    f[k] = y(k) * std::exp(2*M_PI*i*(static_cast<double>(n)/static_cast<double>(N*k))); 
  }
  Eigen::FFT<double> fft;
  fft.fwd(c, f); // -> c = fft(f);

  // From \eqref{eq:trigpcomp}: \Blue{$\alpha_j = \frac{1}{2}(\gamma_{n-j}+\gamma_{n+j})$} and \Blue{$\beta_j = \frac{1}{2i}(\gamma_{n-j}-\gamma_{n+j})$}, \Blue{$j=1,\ldots,n$}, \Blue{$\alpha_0 = \gamma_n$} 
  a = Eigen::VectorXcd(n + 1);
  b = VectorXcd(n);

  a(0) = c[n];
  for (unsigned l = 1; l <= n; ++l) {
    a(l) = c[n - l] + c[n + l];
    b(l - 1) = -i*( c[n - l] - c[n + l] );
  }
  // dont forget scaling factor of forward FFT!
  a /= N; b /= N;
}
/* SAM_LISTING_END_2 */


} //namespace trigipequid
