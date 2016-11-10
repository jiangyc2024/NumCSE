# include <iostream>
# include <cmath> 
# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXd;
using Eigen::VectorXcd; // complex Eigen vector

/* SAM_LISTING_BEGIN_0 */
// Efficient computation of coefficients in expansion \eqref{eq:trigpreal} for a trigonometric
// interpolation polynomial in equidistant points \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}
// IN : \texttt{y} has to be a row vector of odd length, return values are column vectors
//      \texttt{a}, \texttt{b} will be used to save the expansion coefficients
template <class SCALAR>
std::pair<VectorXcd,VectorXcd> trigipequid(const Eigen::Matrix<SCALAR,Eigen::Dynamic,1>& y) {
  const unsigned N = y.size();
  assert (N % 2 == 1 && "Number of points must be odd!");

  const unsigned n = (N - 1)/2;
  // prepare data for fft
  std::complex< double > i(0,1); // imaginary unit
  std::vector< std::complex<double> > f(N);
  std::vector< std::complex<double> > c;

  for (unsigned k = 0; k < N; ++k) {
    // see \eqref{tip:FM}
    f[k] = y(k) * std::exp(2*M_PI*i*(double(n)/N*k)); 
  }
  Eigen::FFT<double> fft;
  fft.fwd(c, f); // -> c = fft(f);

  // From \eqref{eq:trigpcomp}: \Blue{$\alpha_j = \frac{1}{2}(\gamma_{n-j}+\gamma_{n+j})$} and \Blue{$\beta_j = \frac{1}{2i}(\gamma_{n-j}-\gamma_{n+j})$}, \Blue{$j=1,\ldots,n$}, \Blue{$\alpha_0 = \gamma_n$} 
  VectorXcd a = Eigen::VectorXcd(n + 1);
  VectorXcd b = VectorXcd(n);

  a(0) = c[n];
  for (unsigned l = 1; l <= n; ++l) {
    a(l) = c[n - l] + c[n + l];
    b(l - 1) = -i*( c[n - l] - c[n + l] );
  }
  // dont forget scaling factor of forward FFT!
  a /= N; b /= N;

  return std::make_pair(a,b);
}
/* SAM_LISTING_END_0 */
