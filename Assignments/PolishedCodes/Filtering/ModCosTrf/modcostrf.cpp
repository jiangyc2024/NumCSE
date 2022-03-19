/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2021
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <vector>

#define _USE_MATH_DEFINES

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd modcostrf(const Eigen::VectorXd &a) {
  using Comp = std::complex<double>;
  unsigned int n = a.size();
  Eigen::VectorXd f(n);
  Eigen::FFT<double> fft;
  Eigen::VectorXcd tmp(2 * n);
  tmp << 0.5 * a, Eigen::VectorXcd::Zero(n);
  Eigen::VectorXcd v1 = fft.fwd(tmp);
  tmp << 0.5 * a.reverse(), Eigen::VectorXcd::Zero(n);
  Eigen::VectorXcd v2 = fft.fwd(tmp);
  Comp fac = std::exp(((n - 1) * M_PI * Comp(0, 1)) / (double)n);
  Comp d(1.0, 0.0);
  for (int k = 0; k < n; ++k) {
    f[k] = (v1[k] + d * v2[k]).real();
    d *= fac;
  }
  return f;
}
/* SAM_LISTING_END_1 */

Eigen::VectorXd modcostrf_simple(const Eigen::VectorXd &a) {
  unsigned int n = a.size();
  Eigen::VectorXd f(n);
  for (int k = 0; k < n; ++k) {
    f[k] = 0.0;
    for (int j = 0; j < n; ++j) {
      f[k] += a[j] * std::cos((M_PI * j * k) / n);
    }
  }
  return f;
}

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "C++ program for course Numerical Methods for CSE" << std::endl;
  std::cout << "Modified cosine transform via FFT" << std::endl;

  Eigen::VectorXd a = Eigen::VectorXd::LinSpaced(9, 1.0, 9.0);
  std::cout << modcostrf_simple(a).transpose() << std::endl;
  std::cout << modcostrf(a).transpose() << std::endl;

  return 0;
}
