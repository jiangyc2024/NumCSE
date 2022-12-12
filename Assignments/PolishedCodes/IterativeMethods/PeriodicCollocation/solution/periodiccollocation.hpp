#ifndef PERIODICCOLLOCATIONHPP
#define PERIODICCOLLOCATIONHPP
/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>

/**
 * \brief Evaluation of a trigonometric polynomials in M equidistant points
 * $\cob{M\geq N}$
 *
 * \param x coefficient vector
 * \param M number of evaluation points
 * \return Eigen::VectorXd $u_N$
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd eval_uN(const Eigen::VectorXd &x, unsigned int M) {
  const unsigned int N = x.size() - 1;
  assert((M > N) && "Number of evaluation points must be > N");
  Eigen::VectorXd u = Eigen::VectorXd::Zero(M);
  // TODO: (8-16.a) Evaluate uN efficiently, i.e. with complexity O(M log M)
  // START
  // Zero-pad coefficient vector and copy it into a complex-valued vector
  Eigen::VectorXcd xt = Eigen::VectorXcd::Zero(M);
  xt.head(N + 1) = x;
  // Compute DFT of the vector xt
  Eigen::FFT<double> fft;
  Eigen::VectorXcd y = fft.fwd(xt);
  u = y.real();
  // END
  return u;
}
/* SAM_LISTING_END_1 */

/**
 * \brief Efficient evalution of F
 *
 * \param x vector of evaluation points
 * \return Eigen::VectorXd F(x)
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd eval_F(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  Eigen::VectorXd Fx = Eigen::VectorXd::Zero(N + 1);
  // TODO: (8-16.c) Evaluate F at x.
  // START
  Eigen::ArrayXd s = Eigen::ArrayXd::LinSpaced(N + 1, 0, N);
  // Contribution $\cob{\left[\sin(\pi\frac{k}{N+1})\right]_k}$
  const Eigen::VectorXd rhs = ((M_PI / (N + 1)) * s).sin().matrix();
  // Vector $\cob{\left[x_j \cdot 4\pi^2j^2\right]_j}$
  s = 4 * M_PI * M_PI * s * s;
  const Eigen::VectorXd d{(x.array() * s).matrix()};
  Fx = eval_uN(d, N + 1) + (eval_uN(x, N + 1).array().pow(3)).matrix() - rhs;
  // END
  return Fx;
}
/* SAM_LISTING_END_2 */

/**
 * \brief Efficient evaluation of the Jacobian of F
 *
 * \param x vector of evaluation points
 * \return Eigen::MatrixXd Jacobian of F at x
 */
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd eval_DF(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  const Eigen::VectorXd u{eval_uN(x, N + 1)};
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N + 1, N + 1);
  constexpr double fourpisq = 4 * M_PI * M_PI;
  // TODO: (8-16.d) Evaluate the Jacobian of F at x.
  // START
  for (unsigned int k = 0; k <= N; ++k) {
    const double fac = (2 * M_PI * k) / static_cast<double>(N + 1);
    const double uNksq = 3.0 * u[k] * u[k];
    for (unsigned int j = 0; j <= N; ++j) {
      J(k, j) = (fourpisq * j * j + uNksq) * std::cos(fac * j);
    }
  }
  // END
  return J;
}
/* SAM_LISTING_BEGIN_4 */

#endif
