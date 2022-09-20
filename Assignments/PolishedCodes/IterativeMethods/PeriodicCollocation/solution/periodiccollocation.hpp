/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
#ifndef PERIODICCOLLOCATIONHPP
#define PERIODICCOLLOCATIONHPP

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
#include "helper_functions.hpp"

// Evaluation of a trigonometric polynomials in M equidistant points $\cob{M\geq
// N}$
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd eval_uN(const Eigen::VectorXd &x, unsigned int M) {
  unsigned int N = x.size() - 1;
  assert((M > N) && "Number of evaluation points must be > N");
  Eigen::VectorXd u(M);
  // START: Student solution
  // Zero-pad coefficient vector and copy it into a complex-valued vector
  Eigen::VectorXcd xt = Eigen::VectorXcd::Zero(M);
  xt.head(N + 1) = x;
  // Compute DFT of the vector xt
  Eigen::FFT<double> fft;
  Eigen::VectorXcd y = fft.fwd(xt);
  u = y.real();
  // END: student solution
  return u;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd eval_F(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  Eigen::VectorXd Fx;
  // START Student code
  Eigen::ArrayXd s = Eigen::ArrayXd::LinSpaced(N + 1, 0, N);
  // Contribution $\cob{\left[\sin(\pi\frac{k}{N+1})\right]_k}$
  const Eigen::VectorXd rhs = ((M_PI / (N + 1)) * s).sin().matrix();
  // Vector $\cob{\left[x_j \cdot 4\pi^2j^2\right]_j}$
  s = 4 * M_PI * M_PI * s * s;
  const Eigen::VectorXd d{(x.array() * s).matrix()};
  Fx = eval_uN(d, N + 1) + (eval_uN(x, N + 1).array().pow(3)).matrix() - rhs;
  // END Student code
  return Fx;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd eval_DF(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  const Eigen::VectorXd u{eval_uN(x, N + 1)};
  Eigen::MatrixXd J(N + 1, N + 1);
  double fourpisq = 4 * M_PI * M_PI;
  for (unsigned int k = 0; k <= N; ++k) {
    const double fac = (2 * M_PI * k) / static_cast<double>(N + 1);
    const double uNksq = 3.0 * u[k] * u[k];
    for (unsigned int j = 0; j <= N; ++j) {
      J(k, j) = (fourpisq * j * j + uNksq) * std::cos(fac * j);
    }
  }
  return J;
}
/* SAM_LISTING_BEGIN_4 */

#endif
