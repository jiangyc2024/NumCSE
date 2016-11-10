# include <iostream>
# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::ArrayXcd;
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
VectorXcd trigipequidcomp(const VectorXcd& a, const VectorXcd& b, const unsigned N) {
  const unsigned n = a.size() - 1;
  assert (N >= (2*n - 1) && "N is too small! Must be larger than 2*a.size() - 4.");

  const std::complex<double> i(0,1); // imaginary unit
  // build vector \Blue{$\gamma$}
  VectorXcd gamma(2*n + 1);
  gamma(n) = a(0);
  for (unsigned k = 0; k < n; ++k) {
    gamma(k) = 0.5*( a(n - k) + i*b(n - k - 1) );
    gamma(n + k + 1) = 0.5*( a(k + 1) - i*b(k) );
  }

  // zero padding
  VectorXcd ch(N); ch << gamma, VectorXcd::Zero(N - (2*n + 1));

  // build conjugate fourier matrix
  Eigen::FFT<double> fft;
  VectorXcd chCon = ch.conjugate(); 
  VectorXcd v = fft.fwd(chCon).conjugate();

  // multiplicate with conjugate fourier matrix
  VectorXcd y(N);
  for (unsigned k = 0; k < N; ++k) {
    y(k) = v(k) * std::exp( -2.*k*n*M_PI/N*i );
  }

  return y;
}
/* SAM_LISTING_END_0 */
