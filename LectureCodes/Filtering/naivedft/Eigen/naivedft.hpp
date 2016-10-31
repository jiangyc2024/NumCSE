# include <complex>
# include <Eigen/Dense>
using Eigen::VectorXcd;

void naivedft(const VectorXcd& y, VectorXcd& c) {
  const long n = y.size();
  const std::complex<double> i(0,1); 
  c.resize(n);
  std::complex<double> w = std::exp(-2*M_PI/n*i),
                       s = w;

  c(0) = y.sum();
  for (long j = 1; j < n; ++j) {
    c(j) = y(n-1);
    for (long k = n-2; k >= 0; --k) {
      c(j) = c(j)*s + y(k);
    }
    s *= w;
  }
}
