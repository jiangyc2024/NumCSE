/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#include <stdexcept>
#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>

namespace SumRandVar {

/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd pconv(const Eigen::VectorXd &u, const Eigen::VectorXd &x) {
  const int n = x.size();
  assert(u.size() == n);
  Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
  // ``naive'' two-loop implementation of discrete periodic convolution
  for (int k = 0; k < n; ++k) {
    for (int j = 0, l = k; j <= k; ++j, --l) z[k] += u[l] * x[j];
    for (int j = k + 1, l = n - 1; j < n; ++j, --l) z[k] += u[l] * x[j];
  }
  return z;
}
/* SAM_LISTING_END_0 */

// General discrete convolution: double loop implementation
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd dconv(const Eigen::VectorXd &h, const Eigen::VectorXd &x) {
  const int n = h.size();
  const int m = x.size();
  Eigen::VectorXd y{Eigen::VectorXd::Zero(n + m - 1)};
  for (int k = 0; k < n + m - 1; ++k) {
    for (int j = std::max(0, k - n + 1); j < std::min(m, k + 1); ++j) {
      y[k] += h[k - j] * x[j];
    }
  }
  return y;
}
/* SAM_LISTING_END_1 */

// N-fold convolution
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd convpow_loop(const Eigen::VectorXd &p, unsigned int N) {
  // Vector in which multiple convolutions are accumulated
  Eigen::VectorXd acc{p};
  // Carry out N-1 convolutions
  for (unsigned int k = 1; k < N; ++k) {
    const Eigen::VectorXd tmp = dconv(acc, p);
    acc = tmp;
  }
  return acc;
}
/* SAM_LISTING_END_2 */

// N-fold convolution by periodic convolutions
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd convpow_per(const Eigen::VectorXd &p, unsigned int N) {
  const int n = p.size();
  // Zero-padded vector representing M-periodic sequence
  const int M = N * (n - 1) + 1;
  Eigen::VectorXd q_per =
      (Eigen::VectorXd(M) << p, Eigen::VectorXd::Zero(M - n)).finished();
  // Vector in which multiple convolutions are accumulated
  Eigen::VectorXd acc_per{q_per};
  Eigen::VectorXd tmp_per(M);
  // Carry out N-1 convolutions
  for (unsigned int k = 1; k < N; ++k) {
    tmp_per = pconv(acc_per, q_per);
    acc_per = tmp_per;
  }
  return acc_per;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convpow_pdft(const Eigen::VectorXd &p, unsigned int N) {
  const int n = p.size();
  // Zero-padded vector representing M-periodic sequence
  const int M = N * (n - 1) + 1;
  Eigen::VectorXd q_per =
      (Eigen::VectorXd(M) << p, Eigen::VectorXd::Zero(M - n)).finished();
  // Vector in which multiple convolutions are accumulated
  Eigen::VectorXd acc_per{q_per};
  Eigen::VectorXd tmp_per(M);
  // Carry out N-1 convolutions
  for (unsigned int k = 1; k < N; ++k) {
    tmp_per = pconv(acc_per, q_per);
    acc_per = tmp_per;
  }
  return acc_per;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_7 */
Eigen::VectorXd convpow_fft(const Eigen::VectorXd &p, unsigned int N) {
  using Comp = Eigen::VectorXcd::Scalar;
  const unsigned int n = p.size();
  const unsigned int M = N * (n - 1) + 1;
  // Helper object for DFT in Eigen
  Eigen::FFT<double> fft;
  return fft
      .inv(fft.fwd((Eigen::VectorXcd(M) << p.cast<Comp>(),
                    Eigen::VectorXcd::Zero(M - n))
                       .finished())
               .unaryExpr([N](const Comp &z) -> Comp { return std::pow(z, N); })
               .eval())
      .real();
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd convpow(const Eigen::VectorXd &p, unsigned int N) {
  using Comp = Eigen::VectorXcd::Scalar;
  const unsigned int n = p.size();
  // Length of period for periodic embedding
  const unsigned int M = N * (n - 1) + 1;
  // Zero-padded complex vector representing an M-periodic sequence
  Eigen::VectorXcd q_per =
      (Eigen::VectorXcd(M) << p.cast<Comp>(), Eigen::VectorXcd::Zero(M - n))
          .finished();
  // Helper object for DFT in Eigen
  Eigen::FFT<double> fft;
  // Compute discrete Fourier transform
  Eigen::VectorXcd q_hat = fft.fwd(q_per);
  // Take to the N-th power componentwise
  for (unsigned int l = 0; l < q_hat.size(); ++l) {
    q_hat[l] = std::pow(q_hat[l], N);
  }
  // Inverse DFT. Its real part gives the final (N-1)-fold convolution
  return fft.inv(q_hat).real();
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd pconvfft(const Eigen::VectorXd &u, const Eigen::VectorXd &x) {
  const int M = u.size();
  assert(x.size() == M);
  // FFT in eigen requires complex argument vectors
  Eigen::FFT<double> fft;
  const Eigen::VectorXcd u_comp = u.cast<typename Eigen::VectorXcd::Scalar>();
  const Eigen::VectorXcd x_comp = x.cast<typename Eigen::VectorXcd::Scalar>();
  return (fft.inv(((fft.fwd(u_comp)).cwiseProduct(fft.fwd(x_comp))).eval()))
      .real();
}
/* SAM_LISTING_END_4 */

// Probability distribution of the sum of the N random variables
// Very compact implementation using chaining of Eigen operations. 
/* SAM_LISTING_BEGIN_8 */
Eigen::VectorXd probDistSum(const Eigen::VectorXd &p, unsigned int N) {
  using Comp = Eigen::VectorXcd::Scalar;
  if ((p.minCoeff() < 0.0) || (p.maxCoeff() > 1.0) ||
      (std::abs(p.sum() - 1.0) > 1.0E-10)) {
    throw std::runtime_error(
        "p does not conform with a probability distribution");
  }
  const unsigned int n = p.size();
  const unsigned int M = N * (n - 1) + 1;
  Eigen::FFT<double> fft;
  return fft
      .inv(fft.fwd((Eigen::VectorXcd(M) << p.cast<Comp>(),
                    Eigen::VectorXcd::Zero(M - n))
                       .finished())
               .unaryExpr([N](const Comp &z) -> Comp { return std::pow(z, N); })
               .eval())
      .real();
}
/* SAM_LISTING_END_8 */

}  // namespace SumRandVar

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Comvolution by periodic convolution" << std::endl;

  Eigen::VectorXd p(3);
  p << 0.1, 0.5, 0.4;
  std::cout << " p = " << p.transpose() << std::endl;
  std::cout << "conv(p,p) = " << SumRandVar::dconv(p, p).transpose()
            << std::endl;
  std::cout << "pconv(p,p) = " << SumRandVar::pconv(p, p).transpose()
            << std::endl;
  std::cout << "pconvfft(p,p) = " << SumRandVar::pconvfft(p, p).transpose()
            << std::endl;

  for (unsigned int N = 1; N < 6; ++N) {
    std::cout << "convpow: p*^" << N << " = "
              << SumRandVar::convpow(p, N).transpose() << std::endl;
  }

  for (unsigned int N = 1; N < 6; ++N) {
    std::cout << "convpow_per: p*^" << N << " = "
              << SumRandVar::convpow_per(p, N).transpose() << std::endl;
  }

  for (unsigned int N = 1; N < 6; ++N) {
    std::cout << "convpow_pdft: p*^" << N << " = "
              << SumRandVar::convpow_pdft(p, N).transpose() << std::endl;
  }

  for (unsigned int N = 1; N < 6; ++N) {
    std::cout << "convpow_loop: p*^" << N << " = "
              << SumRandVar::convpow_loop(p, N).transpose() << std::endl;
  }

  for (unsigned int N = 1; N < 6; ++N) {
    std::cout << "convpow_fft: p*^" << N << " = "
              << SumRandVar::convpow_fft(p, N).transpose() << std::endl;
  }
  return 0;
}
