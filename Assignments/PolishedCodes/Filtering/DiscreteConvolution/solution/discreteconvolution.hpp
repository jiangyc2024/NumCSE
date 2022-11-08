#ifndef DISCRETECONVOLUTION_HPP
#define DISCRETECONVOLUTION_HPP

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <unsupported/Eigen/FFT>

/**
 * @brief A ``naive'' two-loop implementation of discrete periodic convolution
 *
 * @param u vector of size n
 * @param x vector of size n to be convolved with u
 * @return Eigen::VectorXd discrete periodic convolution of u and x
 */
/* SAM_LISTING_BEGIN_1 */
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
/* SAM_LISTING_END_1 */

/**
 * @brief Efficient discrete periodic convolution
 *
 * @param p vector of size n
 * @param x vector of size n to be convolved with p
 * @return Eigen::VectorXd discrete periodic convolution of p and x
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd pconv_fast(const Eigen::VectorXd &p, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(n == p.size());
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);

  // TODO: (4-5.b) Compute the discrete periodic convolution of p with x
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
/* SAM_LISTING_END_2 */

/**
 * @brief Resized periodic discrete convolution
 *
 * @param p vector of size n
 * @param x vector of size n to be convolved with p
 * @param N
 * @return Eigen::VectorXd discrete periodic convolution of p and x
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd pconvN(const Eigen::VectorXd &p, const Eigen::VectorXd &x,
                       unsigned int N) {
  const int n = x.size();
  assert(n == p.size());
  // Requirement on N to avoid overlap in periodic padding
  assert(N >= 2 * n);
  // Periodic padding of vector p
  Eigen::VectorXd pt{Eigen::VectorXd::Zero(N)};
  pt.head(n) = p;
  pt.tail(n) = p;
  // Zero padding of vector x
  Eigen::VectorXd xt{Eigen::VectorXd::Zero(N)};
  xt.head(n) = x;
  // Periodic convolution of length N
  const Eigen::VectorXd y = pconv(pt, xt);
  return y.head(n);
}
/* SAM_LISTING_END_3 */

/**
 * @brief Discrete convolution by periodic discrete convolution
 *
 * @param h vector of size n
 * @param x vector of size m
 * @return Eigen::VectorXd discrete periodic convolution of h and x
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd dconv_p(const Eigen::VectorXd &h, const Eigen::VectorXd &x) {
  const int n = h.size();   // length of vector $\cob{\Vh}$
  const int m = x.size();   // length of vector $\cob{\Vx}$
  const int N = m + n - 1;  // Minimal length of periodic convolution
  // Zero-padded vectors of length m+n-1
  Eigen::VectorXd ht = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xt = Eigen::VectorXd::Zero(N);
  ht.head(n) = h;
  xt.head(m) = x;
  // Discrete periodic convolution, see \lref{cpp:pconv}
  return pconv(ht, xt);
}
/* SAM_LISTING_END_4 */

#endif