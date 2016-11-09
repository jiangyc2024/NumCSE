# include <iostream>
# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
// Recursive DFT for vectors of length \Blue{$n=2^L$}
VectorXcd fftrec(const VectorXcd& y) {
  using idx_t = VectorXcd::Index;
  using comp = std::complex<double>;
  // Nothing to do for DFT of length 1
  const idx_t n = y.size();
  if (n == 1) return y; 
  if (n % 2 != 0) throw std::runtime_error("size(y) must be even!");
  const Eigen::Map<const Eigen::Matrix<comp,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>
    Y(y.data(),n/2,2); // for selecting even and off compoents
  const VectorXcd c1 = fftrec(Y.col(0)),c2 = fftrec(Y.col(1));
  const comp omega = std::exp(-2*M_PI/n*comp(0,1)); // root of unity \Blue{$\omega_n$}
  comp s(1.0,0.0);
  VectorXcd c(n);
  // Scaling of DFT of odd components plus periodic copying
  for (long k = 0; k < n; ++k) {
    c(k) = c1(k%(n/2)) + c2(k%(n/2))*s;
    s *= omega;
  }
  return c;
}
/* SAM_LISTING_END_0 */

/*
  Eigen::FFT<double> fft;
  VectorXcd c1( n/2 ), c2( n/2 ), yh( n/2 );
  for (long j = 0; j < n/2; ++j)  {
    yh(j) = y(2*j);
    c1 = fft.fwd(yh);
  }
  for (long l = 0; l < n/2; ++l) {
    yh(l) = y(2*l + 1);
    c2 = fft.fwd(yh);
  }
  std::complex<double> i(0,1); // imaginary unit
  c.resize(n);
  for (long k = 0; k < n; ++k) {
    c(k) = c1(k%(n/2)) + c2(k%(n/2))*std::exp(-2*M_PI/n*k*i);
  }
}
*/
