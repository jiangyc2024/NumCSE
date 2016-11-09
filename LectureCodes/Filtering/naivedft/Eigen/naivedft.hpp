# include <complex>
# include <Eigen/Dense>
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
// DFT \eqref{eq:DFT} of vector y returned in c
void naivedft(const VectorXcd& y, VectorXcd& c) {
  using idx_t = VectorXcd::Index;
  const idx_t n = y.size();
  const std::complex<double> i(0,1); 
  c.resize(n);
  // root of unity \Blue{$\omega_n$}, w holds its powers
  std::complex<double> w = std::exp(-2*M_PI/n*i),s = w;
  c(0) = y.sum();
  for (long j = 1; j < n; ++j) {
    c(j) = y(n-1);
    for (long k = n-2; k >= 0; --k) c(j) = c(j)*s + y(k); 
    s *= w;
  }
}
/* SAM_LISTING_END_0 */
