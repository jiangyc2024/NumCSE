
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>

using namespace Eigen;

// For reference, a ``naive'' two-loop implementation
// of discrete periodic convolution.
Eigen::VectorXd pconv(const Eigen::VectorXd &u, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(n == u.size());
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
  for (int k = 0; k < n; ++k) {
    // Calculate z[k] = \sum_{j=0}^{n-1} u[(k-j)%n] * x[j]
    for (int j = 0, l = k; j <= k; ++j, --l) {
      z[k] += u[l] * x[j];
    }
    for (int j = k + 1, l = n - 1; j < n; ++j, --l) {
      z[k] += u[l] * x[j];
    }
  }
  return z;
}

Eigen::VectorXd pconv_fast(const Eigen::VectorXd &p, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(n == p.size());
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);

  // TO DO: (0-2.b) Compute the discrete periodic convolution of p with x
  // in an efficient manner.
  // START

  // N >= 2*n is required to avoid overlap in periodic padding
  int N = 1;
  while (N < 2 * n) {
    N *= 2;
  }
  // Note that Eigen's FFT routines expect complex vectors
  // \cor{Periodic padding} of vector p
  Eigen::VectorXcd pt{Eigen::VectorXcd::Zero(N)};
  pt.head(n) = p.template cast<std::complex<double>>();
  pt.tail(n) = p.template cast<std::complex<double>>();
  // Zero padding of vector x
  Eigen::VectorXcd xt{Eigen::VectorXcd::Zero(N)};
  xt.head(n) = x.template cast<std::complex<double>>();
  // Periodic convolution of length $\cob{N = 2^l}$ realized by FFT
  Eigen::FFT<double> fft;
  z = (fft.inv(((fft.fwd(pt)).cwiseProduct(fft.fwd(xt))).eval()))
      .real()
      .head(n);

  // END
  return z;
}
