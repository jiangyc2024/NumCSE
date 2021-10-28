
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

  // TO DO: Compute the discrete periodic convolution of p with x
  // in an efficient manner.
  // START
 
  // END
  return z;
}
